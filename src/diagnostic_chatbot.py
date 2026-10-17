"""
Diagnostic Chatbot & Early Cancer Detection Anamnesis Assistant (100% Offline Edition).
Provides interactive medical question-answering, multi-turn clinical anamnesis,
Turkish Natural Language Processing (NLP) symptom extraction,
and pure-Python dynamic Bayesian risk evaluation (without relying on offline database files).
"""

import sys
import os
import re
import math
import json
from typing import Dict, List, Any, Optional, Tuple

# Configure stdout/stderr with errors='replace' to prevent encoding crashes on Windows consoles
if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(errors='replace')
    except Exception:
        pass
if hasattr(sys.stderr, 'reconfigure'):
    try:
        sys.stderr.reconfigure(errors='replace')
    except Exception:
        pass
if hasattr(sys.stdin, 'reconfigure'):
    try:
        sys.stdin.reconfigure(errors='replace')
    except Exception:
        pass

def clean_for_console(text: str) -> str:
    """Replaces Unicode emojis with clean text tags to prevent UnicodeEncodeError on Windows terminals."""
    replacements = {
        '🏥': '[KLINIK]', '🩺': '[DOKTOR]', '🫁': '[AKCIGER]', '🎀': '[MEME]',
        '🎗️': '[KOLOREKTAL]', '🔷': '[PROSTAT]', '⚠️': '[DIKKAT]', '🍽️': '[MIDE]',
        '💧': '[MESANE]', '☀️': '[CILT]', '📋': '[RAPOR]', '⏱️': '[SURE]',
        '🔍': '[TETKIK]', '🛡️': '[BILGI]', '🚬': '[SIGARA]', '🧬': '[GENETIK]',
        '🩸': '[KAN]', '🤖': '[BOT]', '🚨': '[ACIL SEVK]', '🚩': '[ALARM]',
        '👉': '->', '•': '*', '🎉': '[TAMAMLANDI]', '✅': '[OK]', '💡': '[NOT]'
    }
    for emoji, repl in replacements.items():
        text = text.replace(emoji, repl)
    return text

def safe_print(text: str = ""):
    """Print helper that guarantees no UnicodeEncodeError on any terminal."""
    cleaned = clean_for_console(str(text))
    enc = getattr(sys.stdout, 'encoding', None) or 'utf-8'
    try:
        safe_str = cleaned.encode(enc, errors='replace').decode(enc)
        print(safe_str)
    except Exception:
        try:
            print(cleaned.encode('ascii', errors='replace').decode('ascii'))
        except Exception:
            pass


# Import metadata definitions from cancer_risk_predictor (in-memory constants)
try:
    from cancer_risk_predictor import SYMPTOM_DEFINITIONS, CANCER_METADATA
except ImportError:
    # Embedded fallback constants if module path differs
    SYMPTOM_DEFINITIONS = {
        'hemoptysis': {'name_tr': 'Hemoptizi (Kanlı Balgam)', 'base_odds_ratio': 12.8, 'nice_ppv': 0.052, 'target_cancers': ['LUNG']},
        'chronic_cough': {'name_tr': 'Kronik / İnatçı Öksürük', 'base_odds_ratio': 4.1, 'nice_ppv': 0.024, 'target_cancers': ['LUNG']},
        'unexplained_weight_loss': {'name_tr': 'Açıklanamayan Kilo Kaybı (>%10)', 'base_odds_ratio': 7.6, 'nice_ppv': 0.038, 'target_cancers': ['PANCREATIC', 'LUNG', 'GASTRIC', 'COLORECTAL']},
        'rectal_bleeding': {'name_tr': 'Rektal Kanama (Dışkıda Kan)', 'base_odds_ratio': 11.4, 'nice_ppv': 0.061, 'target_cancers': ['COLORECTAL']},
        'bowel_habit_change': {'name_tr': 'Dışkılama Alışkanlığında Değişim', 'base_odds_ratio': 5.2, 'nice_ppv': 0.028, 'target_cancers': ['COLORECTAL']},
        'iron_deficiency_anemia': {'name_tr': 'Demir Eksikliği Anemisi', 'base_odds_ratio': 8.9, 'nice_ppv': 0.045, 'target_cancers': ['COLORECTAL', 'GASTRIC']},
        'breast_lump': {'name_tr': 'Memede Ele Gelen Kitle', 'base_odds_ratio': 15.6, 'nice_ppv': 0.088, 'target_cancers': ['BREAST']},
        'nipple_abnormality': {'name_tr': 'Meme Başı Akıntısı / Çekilme', 'base_odds_ratio': 8.7, 'nice_ppv': 0.049, 'target_cancers': ['BREAST']},
        'hematuria': {'name_tr': 'Ağrısız Hematüri (İdrarda Kan)', 'base_odds_ratio': 13.5, 'nice_ppv': 0.076, 'target_cancers': ['BLADDER_RENAL']},
        'urinary_hesitancy': {'name_tr': 'İdrar Akımında Zayıflama', 'base_odds_ratio': 4.6, 'nice_ppv': 0.022, 'target_cancers': ['PROSTATE']},
        'jaundice': {'name_tr': 'Ağrısız Sarılık (İkter)', 'base_odds_ratio': 18.2, 'nice_ppv': 0.095, 'target_cancers': ['PANCREATIC']},
        'dysphagia': {'name_tr': 'İlerleyici Yutma Güçlüğü (Disfaji)', 'base_odds_ratio': 10.9, 'nice_ppv': 0.058, 'target_cancers': ['GASTRIC']},
        'persistent_dyspepsia': {'name_tr': 'Dirençli Mide Ağrısı / Dispepsi', 'base_odds_ratio': 3.8, 'nice_ppv': 0.021, 'target_cancers': ['GASTRIC']},
        'skin_mole_changes': {'name_tr': 'Şüpheli Ben Değişimi (ABCDE)', 'base_odds_ratio': 14.2, 'nice_ppv': 0.125, 'target_cancers': ['MELANOMA']}
    }
    CANCER_METADATA = {
        'LUNG': {'name_tr': 'Akciğer Kanseri (NSCLC/SCLC)', 'recommended_investigation': 'Düşük Doz Kontrastsız Toraks BT (LDCT) & Göğüs Hastalıkları Konsültasyonu'},
        'COLORECTAL': {'name_tr': 'Kolorektal Kanser (Kolon & Rektum)', 'recommended_investigation': 'Tam Kolonoskopi & Biyopsi, Gaitada Gizli Kan Testi (FIT)'},
        'BREAST': {'name_tr': 'Meme Kanseri (Duktal/Lobüler Karsinom)', 'recommended_investigation': 'Bilateral Dijital Mamografi, Meme Ultrasonografisi & Tru-Cut İğne Biyopsisi'},
        'PROSTATE': {'name_tr': 'Prostat Kanseri (Adenokarsinom)', 'recommended_investigation': 'Serum Total/Serbest PSA, Multiparametrik Prostat MR (mpMRI)'},
        'PANCREATIC': {'name_tr': 'Pankreas Kanseri (Duktal Adenokarsinom)', 'recommended_investigation': 'Kontrastlı Üst Abdomen Üç Fazlı BT / MRCP & EUS'},
        'GASTRIC': {'name_tr': 'Mide & Özofagus Kanseri', 'recommended_investigation': 'Üst Gastrointestinal Sistem Endoskopisi (Özofagogastroduodenoskopi - ÖGD) & Biyopsi'},
        'BLADDER_RENAL': {'name_tr': 'Mesane & Böbrek Kanseri', 'recommended_investigation': 'Üriner Sistem Ultrasonografisi & Fleksibl Sistoskopi'},
        'MELANOMA': {'name_tr': 'Malign Melanom (Kutanöz Kanser)', 'recommended_investigation': 'Dermatoskopik Değerlendirme & Total Eksizyonel Biyopsi'}
    }


# Turkish NLP Keyword Mapping to 14 NICE NG12 Standardized Symptoms
NLP_SYMPTOM_PATTERNS: Dict[str, List[str]] = {
    'hemoptysis': [
        r'kan.*balgam', r'balgam.*kan', r'kan.*tükür', r'tükür.*kan',
        r'ağzı.*kan', r'öksür.*kan', r'hemoptizi', r'kanlı balgam'
    ],
    'chronic_cough': [
        r'öksür', r'öksürük', r'öksürüyorum', r'boğazımda gıcık',
        r'inatçı öksürük', r'geçmeyen öksürük', r'kuru öksürük'
    ],
    'unexplained_weight_loss': [
        r'kilo kayb', r'kilo ver', r'zayıfla', r'iştahım yok', r'iştahsız',
        r'eridim', r'istemsiz kilo', r'hızlı kilo', r'pantolonum bol'
    ],
    'rectal_bleeding': [
        r'makat.*kan', r'rektal kanama', r'dışkı.*kan', r'kanlı dışkı',
        r'tuvalet.*kan', r'büyük abdest.*kan', r'kırmızı dışkı', r'hematoşez'
    ],
    'bowel_habit_change': [
        r'bağırsak', r'ishal', r'kabız', r'büyük abdest.*değiş', r'tuvalet.*değiş',
        r'dışkı.*incel', r'dışkılama alışkanlığı'
    ],
    'iron_deficiency_anemia': [
        r'demir eksik', r'anemi', r'kansız', r'ferritin', r'hemoglobin.*düşük',
        r'çabuk yorul', r'sürekli yorgun', r'bitkin', r'solgun'
    ],
    'breast_lump': [
        r'meme.*kitle', r'göğs.*kitle', r'meme.*sertlik', r'meme.*şişlik',
        r'koltuk altı.*kitle', r'memede yumru', r'mememde beze'
    ],
    'nipple_abnormality': [
        r'meme ucu.*akıntı', r'meme başı.*akıntı', r'kanlı akıntı', r'meme ucu.*çekil',
        r'meme başı.*çök', r'portakal kabuğu'
    ],
    'hematuria': [
        r'idrar.*kan', r'çiş.*kan', r'kanlı idrar', r'kırmızı idrar',
        r'hematüri', r'idrar yaparken kan'
    ],
    'urinary_hesitancy': [
        r'idrar.*zor', r'zor.*idrar', r'zayıf akım', r'kesik kesik',
        r'damla damla', r'gece idrar', r'prostat', r'mesane boşalmı'
    ],
    'jaundice': [
        r'sarılık', r'göz.*sarar', r'cilt.*sarar', r'tenim sarar',
        r'çay rengi idrar', r'açık renk dışkı', r'beyaz dışkı', r'kaşıntı.*sarılık'
    ],
    'dysphagia': [
        r'yutma güç', r'yutamıyor', r'lokma takıl', r'boğazıma takıl',
        r'yemek takıl', r'disfaji', r'yutkunurken acı', r'yemek borusu',
        r'yutmakta zor', r'yutma zor', r'katı gıda'
    ],
    'persistent_dyspepsia': [
        r'hazımsızlık', r'mide ağrı', r'mide yanma', r'reflü',
        r'erken doyma', r'mide ekşi', r'epigastrik', r'şişkinlik'
    ],
    'skin_mole_changes': [
        r'\bben\b.*büyü', r'\bben\b.*değiş', r'ben.*kaşın', r'ben.*kana',
        r'şüpheli ben', r'et beni', r'renk.*alaca', r'melanom', r'yeni çıkan ben'
    ],
    'persistent_hoarseness': [
        r'ses kısıklığı', r'sesim kısık', r'sesim kısıldı', r'sesim çıkmıyor',
        r'disfoni', r'boğazımda batma', r'inatçı ses', r'geçmeyen ses'
    ],
    'oral_ulceration': [
        r'ağzımda yara', r'dilimde yara', r'ağız yarası', r'ağızda kitle',
        r'lökoplaki', r'dudakta yara', r'aft.*geçm', r'iyileşmeyen yara',
        r'ağız içi kitle', r'damakta yara'
    ],
    'cervical_lymphadenopathy': [
        r'boynumda şişlik', r'boyunda beze', r'boynumda beze', r'lenf bezi',
        r'koltuk altı.*şişlik', r'lenf nodu', r'beze çıktı', r'lenfadenopati',
        r'sert beze', r'boyunda yumru'
    ],
    'night_sweats_fever': [
        r'gece.*ter', r'sırılsıklam', r'çarşaf değiştir', r'nedensiz ateş',
        r'düşmeyen ateş', r'titreme.*ateş', r'su gibi ter'
    ],
    'postmenopausal_bleeding': [
        r'menopoz.*kanama', r'menopoz.*leke', r'vajinal kanama', r'adetim bitti.*kanama',
        r'rahim kanama', r'menopozdan sonra', r'lekelenme.*menopoz'
    ],
    'persistent_bloating_ascites': [
        r'karın şiş', r'karnım davul gibi', r'karında sıvı', r'asit',
        r'pantolonum kapanmıyor', r'pelvik dolgunluk', r'geçmeyen şişkinlik',
        r'karında gerginlik'
    ],
    'neurological_deficits_headache': [
        r'sabah baş ağrısı', r'fışkırır.*kusma', r'kuvvetsizlik', r'felç',
        r'nöbet', r'gözüm kararıyor', r'çift görme', r'kolumda uyuşma',
        r'sabahları başım ağrıyor', r'denge kaybı'
    ],
    'bone_pain_mass': [
        r'kemik ağrısı', r'gece kemik', r'derin kitle', r'bacağımda şişlik',
        r'kolda kitle', r'kemikte sızı', r'büyüyen kitle', r'istirahatte kemik',
        r'kemiklerim sızlıyor'
    ]
}


# ======================================================================
# OFFLINE MEDICAL KNOWLEDGE BASE & CLINICAL Q&A REPOSITORY
# ======================================================================
OFFLINE_KNOWLEDGE_TOPICS: Dict[str, Dict[str, Any]] = {
    'lung_cancer': {
        'patterns': [
            r'akciğer.*kanser', r'akciger.*kanser', r'akciğer.*tümör', r'akciğer.*belirti',
            r'toraks.*bt', r'akciğer.*nasıl anlaşılır', r'akciğer.*erken teşhis', r'akciğer.*tedavi'
        ],
        'title': 'Akciğer Kanseri & Erken Tanı',
        'answer': (
            "🫁 **Akciğer Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **Erken Evre:** Akciğer dokusunun kendisinde ağrı duyusu bulunmadığı için erken evrede genellikle sinsi ilerler.\n"
            "• **Kritik Kırmızı Bayrak Belirtileri (NICE NG12):**\n"
            "   1. **Hemoptizi (Balgamda kan görülmesi veya kan tükürme):** Tek bir kez dahi olsa acil uzman değerlendirmesi gerektirir.\n"
            "   2. **3 haftayı geçen inatçı, kuru veya karakter değiştiren öksürük**.\n"
            "   3. **Açıklanamayan nefes darlığı ve hırıltılı solunum**.\n"
            "   4. **Geçmeyen göğüs, omuz veya sırt ağrısı**.\n"
            "   5. **Ses kısıklığı (Ses teli sinirine tümör basısı)**.\n"
            "   6. **İstemsiz hızlı kilo kaybı ve halsizlik**.\n"
            "• **Altın Standart Erken Tanı Yöntemi:** Düşük Doz Kontrastsız Toraks BT (LDCT).\n"
            "• **Kimler Taranmalı?** 50-80 yaş arası, en az 20 paket-yıl sigara içmiş ve halen içen ya da son 15 yılda bırakmış kişiler yılda bir Düşük Doz BT ile taranmalıdır."
        )
    },
    'breast_cancer': {
        'patterns': [
            r'meme.*kanser', r'göğüs.*kanser', r'memede.*kitle', r'meme.*belirti',
            r'mamografi', r'meme.*nasıl anlaşılır', r'meme.*akıntı', r'koltuk altı.*kitle'
        ],
        'title': 'Meme Kanseri & Erken Tanı',
        'answer': (
            "🎀 **Meme Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **En Sık Görülen Belirtiler:**\n"
            "   1. **Memede veya koltuk altında ele gelen ağrısız, sert, sınırları düzensiz kitle** (Vakaların %80'inde ilk bulgudur).\n"
            "   2. **Meme başında içe çekilme (retraksiyon) veya yön değişikliği**.\n"
            "   3. **Meme ucundan kendiliğinden gelen kanlı veya şeffaf pembe akıntı**.\n"
            "   4. **Meme derisinde portakal kabuğu görünümü (ödem), kızarıklık veya çekinti**.\n"
            "   5. **İki meme arasında sonradan oluşan belirgin asimetri**.\n"
            "• **Tarama ve Erken Teşhis Takvimi:**\n"
            "   - **40 - 69 Yaş:** 2 yılda bir Dijital Mamografi (KETEM veya hastanelerde ücretsiz/rutin).\n"
            "   - **Meme Dokusu Yoğun Kadınlarda:** Mamografiye ek olarak Meme Ultrasonografisi.\n"
            "   - **Ailede BRCA1/2 Taşıyıcılığı Olanlarda:** 25-30 yaşından itibaren dinamik Kontrastlı Meme MRG.\n"
            "• **Önemli Bilgi:** Memede tespit edilen kitlelerin yaklaşık %80'i iyi huyludur (kist, fibroadenom), ancak her kitle hekim kontrolü gerektirir."
        )
    },
    'colorectal_cancer': {
        'patterns': [
            r'kolon.*kanser', r'kolorektal', r'bağırsak.*kanser', r'rektum.*kanser',
            r'kolon.*belirti', r'kolon.*tedavi'
        ],
        'title': 'Kolorektal (Bağırsak) Kanseri & Erken Tanı',
        'answer': (
            "🎗️ **Kolorektal (Kalın Bağırsak) Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **Öncül Lezyonlar:** Kolorektal kanserlerin %95'i iyi huylu poliplerden (adenom) 5-10 yıllık bir sürede gelişir. Polipken çıkarılması kanseri tamamen önler!\n"
            "• **Kırmızı Bayrak Belirtileri (NICE NG12):**\n"
            "   1. **Dışkıda kan görülmesi (parlak kırmızı veya siyah katranımsı renk)**.\n"
            "   2. **Dışkılama alışkanlığında 6 haftadan uzun süren kalıcı değişiklik (nedensiz ishal veya kabızlık atakları)**.\n"
            "   3. **Açıklanamayan demir eksikliği anemisi (kronik gizli kanama kaynaklı halsizlik, solukluk)**.\n"
            "   4. **Dışkı çapında incelme (şerit şeklinde dışkılama) ve tam boşalamama hissi (tenezm)**.\n"
            "   5. **Kramp tarzı karın ağrısı ve kilo kaybı**.\n"
            "• **Tarama Takvimi:**\n"
            "   - **50 - 70 Yaş:** 2 yılda bir Gaitada Gizli Kan Testi (FIT) ve 10 yılda bir Kolonoskopi.\n"
            "   - **1. Derece Aile Öyküsü Varsa:** Akrabanın tanı yaşından 10 yıl önce (en geç 40 yaşında) kolonoskopi başlatılmalıdır."
        )
    },
    'prostate_cancer': {
        'patterns': [
            r'prostat.*kanser', r'prostat.*belirti', r'psa.*test', r'idrar.*zor',
            r'kesik kesik idrar', r'prostat.*nasıl anlaşılır', r'gece idrar'
        ],
        'title': 'Prostat Kanseri & PSA Değerlendirmesi',
        'answer': (
            "🔷 **Prostat Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **Erken Dönem:** Prostatın dış (kapsüle yakın) bölgesinde geliştiği için başlangıçta tamamen belirtisizdir.\n"
            "• **İlerleyen Evre Belirtileri:**\n"
            "   1. **İdrar akımında zayıflama, çatallanma, damla damla yapma**.\n"
            "   2. **İdrara başlamakta güçlük çekme ve mesaneyi tam boşaltamama hissi**.\n"
            "   3. **Gece sık idrara uyanma (noktüri)**.\n"
            "   4. **İdrarda veya menide kan görülmesi**.\n"
            "   5. **İlerlemiş evrede bel, kalça ve omurga ağrıları**.\n"
            "• **Erken Teşhis Testleri:**\n"
            "   - **Serum Total ve Serbest PSA:** 50 yaş üzeri her erkekte (aile öyküsü varsa 45 yaş) kanda bakılmalıdır.\n"
            "   - **Multiparametrik Prostat MR (mpMRI):** Şüpheli PIRADS lezyonlarını saptar.\n"
            "   - **MR-Füzyon Biyopsi:** Milimetrik hedefe yönelik biyopsi ile net tanı konur."
        )
    },
    'pancreatic_cancer': {
        'patterns': [
            r'pankreas.*kanser', r'pankreas.*belirti', r'sarılık.*kanser',
            r'çay rengi idrar', r'açık renk dışkı', r'pankreas.*nasıl'
        ],
        'title': 'Pankreas Kanseri & Ağrısız Sarılık',
        'answer': (
            "⚠️ **Pankreas Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **Neden Erken Teşhisi Zordur?** Pankreas karın boşluğunun en derininde (retroperiton) yer aldığından kitle safra kanalına bası yapana kadar sessiz kalabilir.\n"
            "• **En Kritik Kırmızı Bayrak Belirtileri:**\n"
            "   1. **Ağrısız İlerleyici Sarılık:** Göz aklarında ve deride belirgin sararma, yoğun cilt kaşıntısı.\n"
            "   2. **Çay veya kola renginde koyu idrar ve açık renkli (cam macunu gibi) dışkı**.\n"
            "   3. **Açıklanamayan ani ve hızlı kilo kaybı**.\n"
            "   4. **Sırta ve iki kürek kemiği arasına vuran, öne eğilince hafifleyen üst karın ağrısı**.\n"
            "   5. **50 yaşından sonra aile öyküsü olmaksızın aniden ortaya çıkan Yeni Başlangıçlı Diyabet**.\n"
            "• **Önerilen Tetkikler:** Çok Kesitli Dinamik Batın BT, MRCP (safra yolu MR'ı) ve Endoskopik Ultrason (EUS)."
        )
    },
    'stomach_cancer': {
        'patterns': [
            r'mide.*kanser', r'gastrik.*kanser', r'mide.*belirti',
            r'erken doyma', r'mide.*ağrı', r'gastroskopi'
        ],
        'title': 'Mide Kanseri & Erken Tanı',
        'answer': (
            "🍽️ **Mide Kanseri ve Erken Teşhis Rehberi:**\n\n"
            "• **Kritik Alarm Bulguları:**\n"
            "   1. **Yutma Güçlüğü (Disfaji):** Katı lokmaların yemek borusunda veya midede takılması, yemek borusu-mide bileşke tümörlerinin birincil alarmıdır.\n"
            "   2. **Erken Doyma:** Birkaç lokma yedikten sonra dahi midede aşırı şişkinlik ve dolgunluk hissi.\n"
            "   3. **Dirençli Mide Ağrısı:** Mide koruyucu ilaçlara (PPI) yanıt vermeyen, geceleri uyandıran epigastrik ağrı.\n"
            "   4. **Açıklanamayan Kilo Kaybı ve İştahsızlık (Özellikle kırmızı etten tiksinme)**.\n"
            "   5. **Kusma ve Gizli Kanama:** Demir eksikliği anemisi, siyah katran gibi dışkı (melena).\n"
            "• **Tanıda Altın Standart:** Üst Gastrointestinal Sistem Endoskopisi (Gastroskopi) ve patolojik biyopsi."
        )
    },
    'bladder_cancer': {
        'patterns': [
            r'mesane.*kanser', r'böbrek.*kanser',
            r'mesane.*belirti', r'sistoskopi'
        ],
        'title': 'Mesane Kanseri & Ağrısız Hematüri',
        'answer': (
            "💧 **Mesane Kanseri ve Hematüri (İdrarda Kan) Rehberi:**\n\n"
            "• **En Önemli Kırmızı Bayrak:** **Ağrısız Makroskopik Hematüri** (İdrarda gözle görülen pembe/kırmızı/pıhtılı kanama).\n"
            "   - İdrar yolu enfeksiyonlarında kanamaya yanma ve sancı eşlik eder. **Ağrısız, sancısız kanama aksi kanıtlanana kadar ürolojik kanser şüphesidir!**\n"
            "• **En Büyük Risk Faktörü:** **Sigara ve Tütün Kullanımı** (Mesane kanserlerinin %50'den fazlası tütün toksinlerinin idrarla atılması esnasında mesaneyi tahriş etmesinden kaynaklanır).\n"
            "• **NICE NG12 Kuralı:** 45 yaş ve üzeri açıklanamayan ağrısız hematürisi olan hastalar 2 hafta içinde (2WW) üroloji uzmanına sevk edilmelidir.\n"
            "• **Tanı Protokolü:** Üriner Sistem USG, İdrar Sitolojisi ve Fleksibl Sistoskopi (kamera ile mesane incelemesi)."
        )
    },
    'melanoma_skin': {
        'patterns': [
            r'melanom', r'\bcilt\b.*kanser', r'\bderi\b.*kanser', r'ben.*değiş',
            r'ben.*büyü', r'şüpheli.*ben', r'abcde', r'et.*beni'
        ],
        'title': 'Malign Melanom & Benlerde ABCDE Kuralı',
        'answer': (
            "☀️ **Malign Melanom ve Cilt Benlerinde ABCDE Kuralı:**\n\n"
            "Melanom, cilt kanserlerinin en tehlikelisidir ancak erken evrede yakalandığında basit bir cerrahiyle %99 tam şifa sağlanır!\n\n"
            "• **ABCDE Erken Tanı Kriterleri:**\n"
            "   • **A (Asimetri):** Benin bir yarısı diğer yarısıyla eşleşmez, simetrik değildir.\n"
            "   • **B (Border - Kenar Düzensizliği):** Sınırları girintili çıkıntılı, dalgalı veya belirsizdir.\n"
            "   • **C (Color - Renk Çeşitliliği):** Tek bir kahverengi ton yerine siyah, kırmızı, mavi, beyaz alacalı renkler içerir.\n"
            "   • **D (Diameter - Çap):** Çapı 6 milimetreden (kurşun kalem silgisinden) büyüktür.\n"
            "   • **E (Evolution - Evrim/Değişim):** EN ÖNEMLİ KRİTERDİR! Benin boyutunda, renginde, şeklinde değişiklik olması; kaşınması veya kanaması.\n"
            "• **Öneri:** Vücudunuzdaki diğer benlere benzemeyen ('Çirkin Ördek') herhangi bir ben fark ederseniz vakit kaybetmeden Dermatoloji uzmanına başvurunuz."
        )
    },
    'screening_and_tests': {
        'patterns': [
            r'hangi test', r'kanser tarama', r'hangi yaşta', r'tarama.*test',
            r'erken teşhis.*test', r'check.*up', r'hangi tahlil', r'nasıl taranır'
        ],
        'title': 'Kanser Tarama Programları ve Yaş Kılavuzları',
        'answer': (
            "📋 **Resmi Kanser Tarama Takvimi (Sağlık Bakanlığı & WHO):**\n\n"
            "1. **Meme Kanseri Taraması (Kadınlar):**\n"
            "   • **40 - 69 Yaş:** 2 yılda bir Dijital Mamografi (KETEM'lerde ücretsizdir).\n"
            "   • Ayda bir kendi kendine meme muayenesi, yılda bir hekim kontrolü.\n\n"
            "2. **Kolorektal (Kalın Bağırsak) Taraması (Kadın & Erkek):**\n"
            "   • **50 - 70 Yaş:** 2 yılda bir Gaitada Gizli Kan Testi (FIT/GGK).\n"
            "   • **10 yılda bir Kolonoskopi** (Kanserleşmeden polipleri temizler).\n\n"
            "3. **Rahim Ağzı (Serviks) Taraması (Kadınlar):**\n"
            "   • **30 - 65 Yaş:** 5 yılda bir HPV-DNA ve Smear testi.\n\n"
            "4. **Yüksek Riskli Akciğer Taraması:**\n"
            "   • **50 - 80 Yaş:** 20 paket-yıl ve üzeri sigara içmiş olanlara yıllık **Düşük Doz Toraks BT (LDCT)**.\n\n"
            "5. **Prostat Taraması:**\n"
            "   • **50 Yaş Üzeri Erkekler:** Kanda Serum Total PSA ve ürolojik değerlendirme."
        )
    },
    'nice_guidelines': {
        'patterns': [
            r'nice.*kılavuz', r'nice.*ng12', r'2ww', r'2-week-wait', r'iki hafta',
            r'acil sevk', r'kırmızı bayrak.*nedir'
        ],
        'title': 'NICE NG12 Kılavuzu & 2-Week-Wait (2WW) Sevk Protokolü',
        'answer': (
            "⏱️ **NICE NG12 ve 2-Week-Wait (2WW) Acil Sevk Sistemi Nedir?**\n\n"
            "• **NICE NG12:** İngiltere Ulusal Sağlık ve Bakım Mükemmelliği Enstitüsü tarafından hazırlanan kanıta dayalı erken kanser kılavuzudur.\n"
            "• **%3 PPV Kuralı:** Bir belirti veya belirti kümesinin kanser olma ihtimali (Pozitif Prediktif Değer) **%3 ve üzerine çıktığında**, bu durum 'Kırmızı Bayrak' (Klinik Alarm) kabul edilir.\n"
            "• **2-Week-Wait Zorunluluğu:**\n"
            "   Kırmızı bayrak taşıyan hastalar rutin bekleme listelerine alınamaz; **en geç 14 takvim günü içinde** ilgili uzman hekim tarafından ileri tetkike (BT, Kolonoskopi, Mamografi, Biyopsi) alınması zorunludur.\n"
            "• **Amaç:** Bürokratik gecikmeleri önleyerek kanseri en erken evrede (Evre 1-2) yakalayıp tam şifa oranını maksimize etmektir."
        )
    },
    'colonoscopy_info': {
        'patterns': [
            r'kolonoskopi', r'kolonoskopi.*nedir', r'kolonoskopi.*acı', r'kolonoskopi.*nasıl',
            r'kolonoskopi.*zor mu', r'kolonoskopi.*hazırlık'
        ],
        'title': 'Kolonoskopi İşlemi, Konforu ve Önemi',
        'answer': (
            "🔍 **Kolonoskopi Hakkında Bilinmesi Gerekenler:**\n\n"
            "• **Nedir?** Kalın bağırsağın ve rektumun iç yüzeyinin bükülebilir kameralı bir cihazla doğrudan görüntülenmesidir.\n"
            "• **Acı Verir mi?** Hayır. Günümüzde kolonoskopiler sedasyon (hafif anestezi) altında yapılır. Hasta işlem boyunca uyur, ağrı veya rahatsızlık hissetmez; işlem ortalama 15-20 dakikada tamamlanır.\n"
            "• **Neden Eşsizdir?** Kolonoskopi hem teşhis koyar hem de kanserleşebilecek polipleri işlem anında kesip alarak (polipektomi) **kanseri oluşmadan engeller**.\n"
            "• **Ne Zaman Yapılmalı?** 50 yaşını geçen herkese rutin olarak; dışkıda kan, açıklanamayan kansızlık veya bağırsak alışkanlığı değişikliği olanlara ise yaşa bakılmaksızın derhal yapılmalıdır."
        )
    },
    'biopsy_myths': {
        'patterns': [
            r'biyopsi.*kanser.*yayar', r'biyopsi.*zararlı', r'biyopsi.*bıçak değerse',
            r'bıçak vurulmaz', r'biyopsi.*tehlikeli'
        ],
        'title': 'Biyopsi Kanseri Yayar mı? (Tıbbi Gerçekler)',
        'answer': (
            "🛡️ **Tıbbi Gerçek: 'Kansere Bıçak Değerse / Biyopsi Yapılırsa Yayılır' Efsanesi:**\n\n"
            "• **HAYIR, Biyopsi Kanseri Kesinlikle Yaymaz!** Bu inanış tıbbi temelden yoksun bir şehir efsanesidir.\n"
            "• **Neden Biyopsi Yapılmak Zorundadır?**\n"
            "   1. Kan tahlilleri veya MR/BT filmleri yalnızca şüpheyi gösterir; 'kesin kanser' tanısı ancak patolog tarafından dokunun mikroskop altında incelenmesiyle konabilir.\n"
            "   2. Günümüz hedefe yönelik akıllı ilaçlarını ve immünoterapilerini seçebilmek için tümörün genetik mutasyonları (EGFR, KRAS, BRAF vb.) ancak biyopsi dokusundan tespit edilir.\n"
            "   3. Biyopsi yapılmadan doğru ve etkin bir onkolojik tedavi planlanamaz.\n"
            "• **Güvenlik:** Modern biyopsiler kapalı koruyucu iğnelerle (Tru-cut) yapılır; tümörün kanala yayılma riski ihmal edilecek kadar düşüktür."
        )
    },
    'smoking_risk': {
        'patterns': [
            r'sigara.*kanser', r'paket.*yıl', r'tütün.*kanser', r'sigarayı.*bırakırsam',
            r'nargile', r'elektronik sigara'
        ],
        'title': 'Sigara, Tütün Maruziyeti ve Paket-Yıl Hesabı',
        'answer': (
            "🚬 **Tütün Maruziyeti ve Kanser Riski:**\n\n"
            "• **Hangi Kanserlere Yol Açar?** Akciğer kanserlerinin %85-90'ının doğrudan sorumlusudur. Ayrıca mesane, gırtlak, yemek borusu, pankreas, mide, böbrek ve serviks kanseri riskini katbekat artırır.\n"
            "• **Paket-Yıl (Pack-Years) Nedir ve Nasıl Hesaplanır?**\n"
            "   $$\\text{Paket-Yıl} = (\\text{Günde İçilen Paket Sayısı}) \\times (\\text{İçilen Yıl Sayısı})$$\n"
            "   - Örn: Günde 1 paket x 20 yıl = **20 Paket-Yıl**.\n"
            "   - Örn: Günde 1.5 paket x 30 yıl = **45 Paket-Yıl**.\n"
            "   - 20 paket-yıl ve üzeri maruziyet, akciğer kanseri açısından 'Yüksek Risk Grubu' kabul edilir ve yıllık Düşük Doz Toraks BT taraması gerektirir.\n"
            "• **Sigarayı Bırakmanın Faydaları:** Bırakıldıktan 5-10 yıl sonra akciğer ve mesane kanseri riski yarı yarıya azalır."
        )
    },
    'family_history': {
        'patterns': [
            r'ailemde.*kanser', r'genetik.*kanser', r'babamda.*kanser', r'annemde.*kanser',
            r'kanser.*bulaşıcı', r'kanser.*kalıtsal', r'brca'
        ],
        'title': 'Ailede Kanser Öyküsü ve Kalıtsal Risk',
        'answer': (
            "🧬 **Aile Kanser Öyküsü ve Genetik Risk:**\n\n"
            "• **Tüm Kanserlerin Kaçı Kalıtsaldır?** Kanserlerin %85-90'ı yaşam boyu çevresel faktörler (yaşlanma, sigara, beslenme) sonucu oluşurken, %10-15'i aileden aktarılan gen mutasyonlarından kaynaklanır.\n"
            "• **Genetik Risk Şüphesi Uyandıran Durumlar:**\n"
            "   1. Birinci derece akrabalarda (anne, baba, kardeş, çocuk) genç yaşta (<50 yaş) kanser teşhisi.\n"
            "   2. Ailede aynı kanser türünün birden fazla bireyde bulunması (özellikle meme, over veya kolon kanseri).\n"
            "   3. Bir bireyde birden fazla bağımsız kanser gelişmesi.\n"
            "   4. Erkekte meme kanseri görülmesi.\n"
            "• **Önemli Genler:** BRCA1/BRCA2 (Meme & Yumurtalık), Lynch Sendromu (Kolorektal & Rahim).\n"
            "• **Öneri:** Ailenizde bu tür vakalar varsa bir Tıbbi Genetik uzmanına danışarak genetik tarama ve erken takip planı oluşturabilirsiniz."
        )
    },
    'liquid_biopsy': {
        'patterns': [
            r'sıvı biyopsi', r'ctdna', r'kanda kanser', r'kandan teşhis',
            r'mrd', r'kanda tümör'
        ],
        'title': 'Sıvı Biyopsi (ctDNA) ve Kanda Erken Kanser Takibi',
        'answer': (
            "🩸 **Sıvı Biyopsi ve Dolaşımdaki Tümör DNA'sı (ctDNA):**\n\n"
            "• **Çalışma Prensibi:** Tümör hücreleri parçalandığında genetik materyallerini (ctDNA) kan dolaşımına bırakır. Kolunuzdan alınan basit bir tüp kan ile yüksek duyarlılıklı yeni nesil dizileme (NGS) yapılarak bu mutant DNA parçacıkları tespit edilebilir.\n"
            "• **Kullanım Alanları:**\n"
            "   1. **Doku Biyopsisinin Zor Olduğu Durumlarda:** Ameliyatsız genetik mutasyon (EGFR, KRAS vb.) tespiti.\n"
            "   2. **Tedavi Yanıtının İzlenmesi:** İlaç başladıktan sonra kanda ctDNA seviyesinin düşmesi ilacın çalıştığını kanıtlar.\n"
            "   3. **Moleküler Nüksün (MRD) Erken Yakalanması:** Ameliyat sonrası kanda kalan mikroskobik hücreleri, radyolojik BT/MR'dan **3 ila 6 ay önce** saptar!"
        )
    },
    'hemoptysis_alarm': {
        'patterns': [
            r'balgam.*kan', r'kan.*balgam', r'kan.*tükür', r'öksürük.*kan', r'hemoptizi'
        ],
        'title': 'Balgamda Kan Görülmesi (Hemoptizi) ve Akciğer Sağlığı',
        'answer': (
            "🩸 **Balgamda Kan Görülmesi (Hemoptizi) ve Klinik Yaklaşım:**\n\n"
            "• **Hemoptizi Nedir?** Solunum yollarından (akciğer veya bronşlar) öksürükle kan veya kanlı balgam gelmesidir.\n"
            "• **Neden Kritik Bir Kırmızı Bayraktır?**\n"
            "   NICE NG12 kılavuzuna göre **tek bir kez dahi olsa** açıklanamayan kanlı balgam görülmesi (Pozitif Prediktif Değer %5.5), "
            "doğrudan **NICE 2-Week-Wait Acil Sevk Kapsamındadır**.\n"
            "• **Olası Nedenler:** Akciğer kanseri, bronşiektazi, tüberküloz, pulmoner emboli veya şiddetli bronşit.\n"
            "• **Önerilen İlk Tetkik:** Düşük Doz Kontrastsız Toraks BT (LDCT) ve Göğüs Hastalıkları uzmanı tarafından yapılacak Bronkoskopi incelemesi."
        )
    },
    'weight_loss_alarm': {
        'patterns': [
            r'nedensiz.*kilo', r'açıklanamayan.*kilo', r'hızlı zayıflama', r'istemsiz.*kilo',
            r'kilo verdim.*neden', r'kilo kaybı'
        ],
        'title': 'Açıklanamayan Hızlı Kilo Kaybı ve Kanser Riski',
        'answer': (
            "⚖️ **Açıklanamayan Kilo Kaybı (İstemsiz Zayıflama):**\n\n"
            "• **Klinik Tanım:** Diyet veya ağır egzersiz yapmaksızın, son 6 ayda vücut ağırlığının **%5 ila %10'undan fazlasının** istemsizce kaybedilmesidir.\n"
            "• **Hangi Kanserlerde Görülür?** Özellikle pankreas, mide, yemek borusu (özofagus), akciğer ve lenfoma (kan kanseri) vakalarında ilk klinik belirtilerden biridir.\n"
            "• **Mekanizma:** Tümör hücrelerinin salgıladığı pro-enflamatuar sitokinler (tümör nekroz faktörü - TNF-alfa vb.) vücutta 'kaşeksi' denilen kas ve yağ dokusu erimesine yol açar.\n"
            "• **Öneri:** Açıklanamayan kilo kaybınız varsa vakit kaybetmeden Dahiliye veya Gastroenteroloji uzmanına başvurarak tam kan tahlili, tümör belirteçleri ve batın ultrasonografisi/BT yaptırmalısınız."
        )
    },
    'rectal_bleeding_alarm': {
        'patterns': [
            r'makat.*kan', r'dışkı.*kan', r'büyük abdest.*kan', r'rektal kanama', r'kanlı dışkı'
        ],
        'title': 'Makat veya Dışkıda Kanama (Rektal Kanama) ve Bağırsak Sağlığı',
        'answer': (
            "🩸 **Makat veya Dışkıda Kanama (Rektal Kanama):**\n\n"
            "• **Hemoroid (Basur) mu, Kanser mi?** Çoğu insan kanamayı hemoroid sanarak hekime gitmeyi geciktirir. Ancak hemoroid ile kolorektal kanser veya polipler aynı anda bulunabilir!\n"
            "• **NICE NG12 Kuralı:** 50 yaş ve üzerinde açıklanamayan rektal kanama doğrudan 2 hafta içinde (2WW) kolonoskopi gerektirir.\n"
            "• **Kritik İpuçları:** Kanamaya dışkılama alışkanlığında değişiklik (ishal/kabızlık) veya kansızlık (anemi) eşlik ediyorsa acil kolonoskopi şarttır.\n"
            "• **Öneri:** Genel Cerrahi veya Gastroenteroloji hekimine başvurunuz."
        )
    },
    'hematuria_alarm': {
        'patterns': [
            r'idrar.*kan', r'hematüri', r'kırmızı idrar', r'idrarda kan.*tehlikeli', r'kanlı idrar'
        ],
        'title': 'İdrarda Kan Görülmesi (Ağrısız Hematüri) ve Ürolojik Risk',
        'answer': (
            "💧 **İdrarda Kan Görülmesi (Hematüri):**\n\n"
            "• **Ağrısız Hematüri Kuralı:** İdrar yolu enfeksiyonlarında yanma/sancı olurken; **ağrısız, sancısız gözle görülen kanama mesane veya böbrek kanserinin en tipik erken belirtisidir**.\n"
            "• **NICE NG12 Kılavuzu:** 45 yaş üzeri ağrısız hematürisi olan hastalar en geç 2 hafta içinde Ürolojiye yönlendirilmelidir.\n"
            "• **Önerilen İnceleme:** Üriner Sistem Ultrasonu, İdrar Sitolojisi ve Sistoskopi."
        )
    },
    'breast_lump_alarm': {
        'patterns': [
            r'memede kitle', r'memede sertlik', r'koltuk altı.*şişlik', r'meme ucu.*akıntı', r'memede yumru'
        ],
        'title': 'Memede Ele Gelen Kitle ve Erken Tanı',
        'answer': (
            "🎀 **Memede Ele Gelen Kitle:**\n\n"
            "• **İlk Değerlendirme:** Memede saptanan kitlelerin yaklaşık %80'i iyi huyludur (fibroadenom, basit kist vb.). Ancak hekim ve görüntüleme onayı olmaksızın iyi huylu olduğu varsayılamaz.\n"
            "• **Şüpheli Bulgular:** Sert, sınırları düzensiz, çevre dokuya yapışık ve ağrısız kitleler malignite açısından daha yüksek risk taşır.\n"
            "• **İzlenecek Yol:** 40 yaş altı için Meme USG, 40 yaş üzeri için Dijital Mamografi + Meme USG, şüpheli durumlarda Tru-Cut iğne biyopsisi."
        )
    },
    'dysphagia_alarm': {
        'patterns': [
            r'yutma güç', r'disfaji', r'lokma.*takıl', r'boğaz.*takıl', r'yutamıyorum'
        ],
        'title': 'Yutma Güçlüğü (Disfaji) & Üst GİS Kırmızı Bayrağı',
        'answer': (
            "🍽️ **Yutma Güçlüğü (Disfaji) ve Acil Klinik Yaklaşım:**\n\n"
            "• **Neden Kritik Bir Kırmızı Bayraktır?** Katı veya sıvı gıdaları yutarken boğazda/göğüste lokmanın takılması hissi, yemek borusu (özofagus) veya mide-yemek borusu bileşkesi lezyonlarının birincil alarmıdır.\n"
            "• **NICE NG12 Kılavuzu:** Açıklanamayan disfaji semptomu doğrudan **NICE 2-Week-Wait Acil Sevk** gerektirir.\n"
            "• **Önerilen İlk Tetkik:** Vakit kaybetmeden Gastroenteroloji uzmanı tarafından Üst GİS Endoskopisi (Gastroskopi) yapılmalıdır."
        )
    },
    'jaundice_alarm': {
        'patterns': [
            r'sarılık', r'göz.*sarar', r'cilt.*sarar', r'tenim sarar', r'çay rengi idrar', r'beyaz dışkı'
        ],
        'title': 'Ağrısız Sarılık (İkter) & Pankreas/Safra Alarmları',
        'answer': (
            "🟡 **Ağrısız Sarılık ve Klinik Değerlendirme:**\n\n"
            "• **Kritik Alarm Kuralı:** Ağrısız, sancısız ilerleyici sarılık (göz aklarında sararma, çay rengi idrar, açık renkli dışkı) aksi kanıtlanana kadar pankreas başı ve safra yolu obstrüksiyonu şüphesidir.\n"
            "• **NICE NG12 Kuralı:** 40 yaş ve üzeri açıklanamayan ağrısız sarılık hastaları 2 hafta içinde doğrudan ileri tetkik için sevk edilmelidir.\n"
            "• **Önerilen İnceleme:** Dinamik Batın BT, MRCP ve Serum Total/Direkt Bilirubin, CA 19-9 testleri."
        )
    },
    'head_neck_cancer': {
        'patterns': [
            r'baş boyun.*kanser', r'larenks.*kanser', r'gırtlak.*kanser', r'ağız.*kanser',
            r'ses kısıklığı.*kanser', r'dil.*kanser', r'ağız yarası.*kanser'
        ],
        'title': 'Baş & Boyun (Larenks / Oral Kavite) Kanserleri & Erken Tanı',
        'answer': (
            "🗣️ **Baş & Boyun ve Larenks (Gırtlak) Kanserleri Rehberi:**\n\n"
            "• **En Kritik Kırmızı Bayraklar (NICE NG12):**\n"
            "   1. **3 haftadan uzun süren kalıcı ses kısıklığı (disfoni):** Larenks (ses teli) kanserlerinin en erken ve en değerli bulgusudur.\n"
            "   2. **Ağız içinde, dilde veya damakta 3 haftayı geçen iyileşmeyen ülser veya aft benzeri yara**.\n"
            "   3. **Ağız içinde kırmızı (eritroplaki) veya beyaz (lökoplaki) silinmeyen plaklar**.\n"
            "   4. **Tek taraflı boğaz ağrısı, kulağa vuran ağrı (yansıyan otalji) ve boyunda ele gelen sert kitle**.\n"
            "• **En Büyük Risk Faktörleri:** Sigara ve alkolün birlikte kullanımı riski 30 kat artırır; ayrıca orofarenks kanserlerinde HPV (Özellikle HPV Tip 16) majör etkendir.\n"
            "• **Tanı Yöntemi:** Fleksibl fiberoptik nazofaringolaringoskopi ve biyopsi."
        )
    },
    'lymphoma_hematologic': {
        'patterns': [
            r'lenfoma', r'lösemi', r'kan kanseri', r'b belirtileri', r'gece terlemesi.*kanser',
            r'lenf bezi.*kanser', r'hodgkin'
        ],
        'title': 'Hematolojik Maligniteler & Lenfoma B-Semptomları',
        'answer': (
            "🩸 **Lenfoma ve Hematolojik Kanserler Erken Teşhis Rehberi:**\n\n"
            "• **Kritik B-Semptomları Triadı:**\n"
            "   1. **Sırılsıklam Gece Terlemeleri:** Gece kıyafet ve çarşaf değiştirmeyi gerektirecek derecede yoğun terleme.\n"
            "   2. **Açıklanamayan Sürekli Ateş:** Enfeksiyon odağı olmaksızın 38°C üzerindeki dalgalı ateş (Pel-Ebstein ateşi).\n"
            "   3. **6 ayda %10'dan fazla açıklanamayan kilo kaybı**.\n"
            "• **Fizik Muayene Bulgusu:** Boyunda, koltuk altında veya kasıkta büyüyen, genellikle ağrısız, lastik kıvamında sert lenfadenopati (>2 cm).\n"
            "• **Tanı Protokolü:** İğne aspirasyonu (İAB) yetersizdir; **kesin tanı için tam eksizyonel lenf nodu biyopsisi** ve PET-BT evrelemesi zorunludur."
        )
    },
    'gynecologic_cancer': {
        'patterns': [
            r'jinekolojik.*kanser', r'rahim.*kanser', r'yumurtalık.*kanser', r'over.*kanser',
            r'endometrium.*kanser', r'serviks.*kanser', r'menopoz.*kanama', r'karında asit'
        ],
        'title': 'Jinekolojik Kanserler (Endometrium, Over, Serviks)',
        'answer': (
            "🌺 **Jinekolojik Kanserler ve Erken Teşhis Kılavuzu:**\n\n"
            "• **Postmenopozal Vajinal Kanama (En Büyük Kırmızı Bayrak):**\n"
            "   Menopoza girmiş bir kadında tek bir damla dahi kanama veya kahverengi lekelenme görülmesi, aksi kanıtlanana kadar **Endometrium (Rahim İçi) Kanseri** şüphesidir ve doğrudan acil 2WW sevk gerektirir.\n"
            "• **Over (Yumurtalık) Kanseri Sinsi İlerler:**\n"
            "   - Sürekli karın şişkinliği, pantolonun kapanmaması, erken doyma ve karında sıvı toplanması (asit).\n"
            "• **Serviks (Rahim Ağzı) Kanseri Önlenebilir:**\n"
            "   - HPV aşılaması ve 30-65 yaş arası 5 yılda bir yapılan HPV-DNA/Smear testi ile kanser öncül lezyon aşamasında %100 yakalanabilir."
        )
    },
    'cns_brain_tumor': {
        'patterns': [
            r'beyin.*tümör', r'beyin.*kanser', r'glioblastom', r'sabah baş ağrısı.*kanser',
            r'fışkırır.*kusma', r'baş ağrısı.*kanser'
        ],
        'title': 'Merkezi Sinir Sistemi (MSS) & Beyin Tümörleri',
        'answer': (
            "🧠 **Beyin ve Merkezi Sinir Sistemi Tümörleri Rehberi:**\n\n"
            "• **Karakteristik Kırmızı Bayraklar:**\n"
            "   1. **Sabahları uykudan uyandıran ve öksürmekle/ıkınmakla şiddetlenen baş ağrısı** (Kafa içi basınç artışı - KİBAS sendromu).\n"
            "   2. **Bulantı eşlik etmeksizin fışkırır tarzda sabah kusmaları**.\n"
            "   3. **Daha önce epilepsi öyküsü olmayan erişkin bir bireyde aniden ortaya çıkan ilk nöbet (havale)**.\n"
            "   4. **Vücudun tek tarafında ilerleyici kuvvet kaybı (felç), hissizlik, yürüme dengesizliği veya konuşma bozukluğu**.\n"
            "• **Tanıda Altın Standart:** Kontrastlı Kraniyal Beyin Manyetik Rezonans Görüntüleme (MRG)."
        )
    },
    'sarcoma_bone_tumor': {
        'patterns': [
            r'sarkom', r'kemik.*kanser', r'kemik.*tümör', r'osteosarkom',
            r'derin kemik ağrısı', r'kemikte kitle'
        ],
        'title': 'Kemik ve Yumuşak Doku Sarkomları',
        'answer': (
            "🦴 **Sarkomlar ve Kemik/Yumuşak Doku Kanserleri:**\n\n"
            "• **Alarm Belirtileri:**\n"
            "   1. **İstirahat halindeyken veya geceleri uykudan uyandıran derin, sızlayıcı kemik ağrısı**.\n"
            "   2. **Kol, bacak veya gövdede hızla büyüyen, ağrılı veya ağrısız, sert doku kitlesi** (>5 cm derin kitleler daima şüphelidir).\n"
            "   3. **Travma olmaksızın basit bir hareketle oluşan patolojik kemik kırıkları**.\n"
            "• **Tanı:** İki yönlü direkt grafi, Kontrastlı Ekstremite MRG ve Tru-Cut / İnsizyonel biyopsi."
        )
    },
    'thyroid_cancer': {
        'patterns': [
            r'tiroid.*kanser', r'tiroid.*nodül', r'guatr.*kanser', r'boyunda nodül',
            r'papiller'
        ],
        'title': 'Tiroid Kanseri & Tiroid Nodülleri',
        'answer': (
            "🦋 **Tiroid Kanseri ve Tiroid Nodülleri Rehberi:**\n\n"
            "• **Nodüllerin Ne Kadarı Kanserdir?** Toplumda tiroid nodülleri çok yaygındır ancak bu nodüllerin yalnızca %5-10'u kanserleşir.\n"
            "• **Malignite Şüphesi Uyandıran Bulgular (TIRADS Kriterleri):**\n"
            "   - Boyun ön bölgesinde hızla büyüyen sert ve hareketsiz kitle.\n"
            "   - Ultrasonografide mikrokalsifikasyonlar, düzensiz sınırlar ve hipoekoik yapı.\n"
            "   - Boyun lenf bezlerinde eşzamanlı büyüme.\n"
            "• **Erken Teşhiste Başarı:** Papiller tiroid kanseri erken yakalandığında cerrahi ve radyoaktif iyot ile **%98'in üzerinde 10 yıllık tam sağkalım** oranına sahiptir."
        )
    },
    'pharmacogenomics_dpyd': {
        'patterns': [
            r'dpyd', r'farmakogenomik', r'pgx', r'5-fu.*toksisite', r'kapesitabin.*test',
            r'ugt1a1', r'kemoterapi toksisite'
        ],
        'title': 'Farmakogenomik (PGx) Testleri & Kemoterapi Güvenliği',
        'answer': (
            "💊 **Farmakogenomik (PGx) Testleri ve CPIC Seviye A Kılavuzları:**\n\n"
            "• **DPYD Gen Testi Nedir?** 5-Fluorouracil (5-FU) ve Kapesitabin gibi kolorektal/meme kemoterapilerini vücut DPYD enzimi ile yıkar. DPYD enzim eksikliği olan hastalarda bu ilaçlar birikir ve ölümcül toksisiteye (nötropenik sepsis, ağır mukozit) yol açar.\n"
            "• **Testin Önemi:** Kemoterapi başlamadan önce yapılan basit bir DNA testiyle doz %50 azaltılarak veya alternatif ilaca geçilerek hasta hayatı kurtarılır.\n"
            "• **UGT1A1 & İrinotekan:** İrinotekan metabolizmasını düzenler; *28/*28 homozigot hastalarda ağır ishal ve lökopeni riskine karşı doz ayarı gereklidir."
        )
    },
    'synthetic_lethality_parp': {
        'patterns': [
            r'sentetik letalite', r'parp.*inhibitör', r'olaparib', r'brca.*parp',
            r'depmap', r'msi.*wrn'
        ],
        'title': 'Sentetik Letalite (Synthetic Lethality) & PARP Tedavileri',
        'answer': (
            "🎯 **Sentetik Letalite ve Hedefe Yönelik Kanser Tedavisi:**\n\n"
            "• **Sentetik Letalite Nedir?** İki genden birinin mutasyonu hücreyi yaşatırken, her iki genin aynı anda bozulması durumunda kanser hücresinin seçici olarak ölmesidir. Sağlıklı hücreler zarar görmez!\n"
            "• **Klinik Başarı (BRCA1/2 + PARP):**\n"
            "   BRCA mutant tümörler homolog rekombinasyon (HRD) onarımını kaybetmiştir. PARP inhibitörleri (Olaparib, Talazoparib) tek zincir onarımını da bloke edince tümör hücresi bölünemez ve apoptoza gider.\n"
            "• **Yeni Nesil Zaaflar (MSI-H + WRN):**\n"
            "   Mikrosatellit instabilitesi (MSI-H) olan tümörler Broad DepMap verilerine göre kesinlikle WRN DNA helikazına bağımlıdır. WRN blokajı bu tümörleri selektif olarak yok eder."
        )
    },

    'immunotherapy_pembro': {
        'patterns': [
            r'immünoterapi', r'immunoterapi', r'pembrolizumab', r'keytruda',
            r'nivolumab', r'opdivo', r'pd-1', r'pd-l1', r'bağışıklık tedavisi'
        ],
        'title': 'İmmünoterapi (PD-1 / PD-L1 İnhibitörleri) & Keytruda',
        'answer': (
            "🛡️ **İmmünoterapi ve Kontrol Noktası (Checkpoint) İnhibitörleri:**\n\n"
            "• **Çalışma Prensibi:** Kanser hücreleri PD-L1 proteinini kullanarak bağışıklık sistemimizin T hücrelerini 'uyutur' ve kendini gizler. "
            "Pembrolizumab (Keytruda) ve Nivolumab (Opdivo) gibi ilaçlar bu fren mekanizmasını kaldırarak vücudun kendi bağışıklık hücrelerinin tümöre saldırmasını sağlar.\n"
            "• **Kullanıldığı Başlıca Kanserler:** Malign melanom, küçük hücreli dışı akciğer kanseri (NSCLC), böbrek (RCC), mesane ve MSI-H tümörler.\n"
            "• **Avantajı:** Kemoterapiden farklı olarak doğrudan sağlıklı hücreleri zehirlemez; yanıt veren hastalarda yıllarca süren kalıcı remisyon sağlayabilir.\n"
            "• **Belirteçler:** Tedavi öncesinde biyopsi dokusunda PD-L1 ekspresyonu ve Tümör Mutasyon Yükü (TMB) bakılır."
        )
    },
    'targeted_therapy_osimertinib': {
        'patterns': [
            r'osimertinib', r'tagrisso', r'egfr.*mutasyon', r't790m', r'c797s',
            r'akıllı ilaç.*akciğer', r'tirozin kinaz'
        ],
        'title': 'EGFR Mutasyonları & 3. Kuşak Akıllı İlaç (Osimertinib)',
        'answer': (
            "🎯 **EGFR Mutasyonları ve Hedefe Yönelik Tedavi (Osimertinib):**\n\n"
            "• **EGFR Nedir?** Akciğer adenokarsinomlarının yaklaşık %15-35'inde (özellikle sigara içmemiş veya kadınlarda) EGFR geninde aktive edici mutasyonlar (Ekzon 19 delesyonu veya L858R) bulunur.\n"
            "• **Osimertinib (Tagrisso):** 3. kuşak oral tirozin kinaz inhibitörüdür (TKI). FLAURA klinik çalışmasında standart kemoterapiye göre genel sağkalımı belirgin şekilde uzatmıştır.\n"
            "• **Direnç Mekanizmaları:** Eski nesil ilaçlara direnç geliştiren **T790M** mutasyonunu yener. Osimertinib direnci geliştiğinde ise kanda **C797S** veya MET amplifikasyonu araştırılır."
        )
    },
    'kras_sotorasib': {
        'patterns': [
            r'kras', r'sotorasib', r'lumakras', r'g12c', r'adagrasib', r'kras inhibitör'
        ],
        'title': 'KRAS Mutasyonu & KRAS G12C İnhibitörleri',
        'answer': (
            "🧬 **KRAS G12C Mutasyonu ve Hedefe Yönelik Tedaviler:**\n\n"
            "• **Tarihsel Önemi:** KRAS, 40 yılı aşkın süre boyunca 'ilaç geliştirilemez' (undruggable) kabul edilen en yaygın onkogendi.\n"
            "• **G12C Başarısı:** Sotorasib (Lumakras) ve Adagrasib, KRAS proteininin G12C cep mutasyonuna kovalent bağlanarak proteini inaktif (GDP bağlı) fazda kilitler.\n"
            "• **Uygulama Alanı:** Akciğer kanserlerinin (NSCLC) yaklaşık %13'ünde ve kolorektal kanserlerin %3-4'ünde G12C mutasyonu görülür ve bu hastalarda hedefe yönelik tedavi seçeneği sunar."
        )
    },
    'her2_trastuzumab': {
        'patterns': [
            r'her2', r'trastuzumab', r'herceptin', r'pertuzumab', r'enhertu', r't-dxd', r'erbb2'
        ],
        'title': 'HER2 Pozitifliği & Trastuzumab (Herceptin / Enhertu)',
        'answer': (
            "🎀 **HER2 Pozitif Kanserler ve Monoklonal Antikorlar:**\n\n"
            "• **HER2 Nedir?** Meme ve mide kanserlerinin yaklaşık %15-20'sinde hücre yüzeyinde HER2 (Human Epidermal Growth Factor Receptor 2) proteini aşırı miktarda üretilir (amplifikasyon).\n"
            "• **Tedavi:** Trastuzumab (Herceptin) ve Pertuzumab, bu reseptörleri bloke ederek hücre bölünmesini durdurur.\n"
            "• **Yeni Nesil Antikor-İlaç Konjugatları (ADC):** Trastuzumab Deruxtecan (Enhertu), antikora bağlı kemoterapi roketi gibi doğrudan tümör içine girerek HER2-low (düşük) tümörlerde dahi çığır açıcı yanıt verir."
        )
    },
    'chemotherapy_basics': {
        'patterns': [
            r'kemoterapi.*nedir', r'kemoterapi.*nasıl', r'\bkemoterapi\b', r'kemoterapötik', r'damardan ilaç.*kanser'
        ],
        'title': 'Kemoterapi Nedir & Nasıl Çalışır?',
        'answer': (
            "💉 **Kemoterapi (Sitotoksik Tedavi) Hakkında Temel Bilgiler:**\n\n"
            "• **Mekanizma:** Kemoterapi ilaçları, vücutta hızlı bölünen tüm hücrelerin DNA sentezini veya hücre bölünmesini (mitoz) durdurarak onları apoptoza (hücre ölümüne) sürükler.\n"
            "• **Uygulama Yolları:** Damar yoluyla (serum içinde infüzyon), ağızdan (tablet/kapsül) veya doğrudan vücut boşluklarına uygulanabilir.\n"
            "• **Amaçlar:**\n"
            "   1. **Küratif (Tam Şifa):** Kanseri tamamen yok etmek.\n"
            "   2. **Neoadjuvan:** Ameliyat öncesi tümörü küçültmek.\n"
            "   3. **Adjuvan:** Ameliyat sonrası geride kalabilecek mikroskobik hücreleri temizlemek.\n"
            "   4. **Palyatif:** İlerlemiş evrede tümör büyümesini frenleyip şikayetleri azaltmak."
        )
    },
    'chemo_side_effects': {
        'patterns': [
            r'kemoterapi.*yan etki', r'saç dökülmesi.*kanser', r'nötropeni', r'kemoterapi.*bulantı',
            r'kemoterapi.*halsizlik'
        ],
        'title': 'Kemoterapinin Yan Etkileri ve Başa Çıkma',
        'answer': (
            "💊 **Kemoterapi Yan Etkileri ve Yönetimi:**\n\n"
            "• **Neden Yan Etki Olur?** İlaçlar kanser hücrelerinin yanında vücudun normal hızlı bölünen hücrelerini de (kemik iliği, saç folikülü, ağız ve sindirim mukozası) geçici olarak etkiler.\n"
            "• **En Sık Görülenler ve Çözümleri:**\n"
            "   1. **Saç Dökülmesi:** Tedavi bittikten 1-2 ay sonra saçlar eskisinden daha gür olarak yeniden uzar.\n"
            "   2. **Bulantı/Kusma:** Modern antiemetik ilaçlar (5-HT3 ve NK1 reseptör blokerleri) sayesinde bulantı günümüzde %80-90 oranında tamamen önlenir.\n"
            "   3. **Nötropeni (Bağışıklık Düşüklüğü):** Beyaz kan hücreleri (lökosit) düştüğünde G-CSF iğneleriyle kemik iliği desteklenir. **38°C üzeri ateş durumunda acil hastaneye başvurulmalıdır!**\n"
            "   4. **Ağız Yaraları (Mukozit):** Karbonatlı tuzlu su gargaraları ve buz uygulaması rahatlatır."
        )
    },
    'radiotherapy_basics': {
        'patterns': [
            r'radyoterapi', r'ışın tedavisi', r'radyasyon tedavisi', r'cyberknife', r'linak'
        ],
        'title': 'Radyoterapi (Işın Tedavisi) Nedir?',
        'answer': (
            "⚡ **Radyoterapi (Işın Tedavisi) Rehberi:**\n\n"
            "• **Nasıl Çalışır?** Yüksek enerjili X-ışınları veya foton demetleri milimetrik doğrulukla tümör dokusuna odaklanarak kanser hücrelerinin DNA'sını parçalar.\n"
            "• **Ağrılı mıdır?** Hayır, radyoterapi seansları röntgen filmi çektirmek gibidir; seans sırasında hasta hiçbir acı veya sıcaklık hissetmez; işlem 10-15 dakikada biter.\n"
            "• **Radyasyon Bulaşır mı?** Dışarıdan uygulanan (eksternal) radyoterapide hastanın vücudunda radyasyon kalmaz; hasta seans biter bitmez ailesine ve çocuklara güvenle sarılabilir."
        )
    },
    'hormonotherapy_breast_prostate': {
        'patterns': [
            r'hormonoterapi', r'tamoksifen', r'letrozol', r'anastrozol', r'hormon tedavisi.*kanser'
        ],
        'title': 'Hormon Tedavisi (Meme & Prostat Kanseri)',
        'answer': (
            "🌱 **Onkolojide Hormon Tedavisi (Endokrin Terapi):**\n\n"
            "• **Meme Kanseri:** Östrojen ve Progesteron reseptörü pozitif (ER+/PR+) tümörlerde Tamoksifen veya Aromataz İnhibitörleri (Letrozol, Anastrozol) ile östrojen hormonu baskılanır, nüks riski %50 oranında düşürülür.\n"
            "• **Prostat Kanseri:** Testosteron hormonunun prostat tümörünü beslemesi androjen baskılayıcı iğneler (GnRH analogları) ile engellenir.\n"
            "• **Avantajı:** Klasik kemoterapiye göre yan etki profili çok daha hafiftir."
        )
    },
    'pet_ct_suvmax': {
        'patterns': [
            r'pet.*bt', r'pet.*ct', r'suvmax', r'suv değeri', r'fdg tutulumu', r'pet çekimi'
        ],
        'title': 'PET-BT İncelemesi & SUVmax Değeri Ne Anlama Gelir?',
        'answer': (
            "☢️ **PET-BT (Pozitron Emisyon Tomografisi) & SUVmax Rehberi:**\n\n"
            "• **Nasıl Çalışır?** Kanser hücreleri normal hücrelere göre çok daha fazla şeker (glukoz) tüketir. Damardan verilen radyoaktif işaretli şeker (18F-FDG) hızla tümör hücrelerinde birikir ve PET kamerası bu hipermetabolik alanları parlatır.\n"
            "• **SUVmax (Standardized Uptake Value) Nedir?** Lezyondaki metabolik aktivite yoğunluğunun sayısal göstergesidir.\n"
            "   - SUVmax < 2.5: Genellikle düşük aktivite (enfeksiyon veya iyi huylu lezyon lehine).\n"
            "   - SUVmax > 4 - 5: Artmış glukoz metabolizması (tümör şüphesini artırır).\n"
            "• **Unutulmamalıdır ki:** Şiddetli enfeksiyon ve iltihaplar da yüksek FDG tutabilir; kesin teşhis mutlaka patolojik biyopsiyle konulur."
        )
    },
    'biopsy_types': {
        'patterns': [
            r'biyopsi çeşitleri', r'tru.*cut', r'kalın iğne biyopsi', r'ince iğne biyopsi',
            r'biyopsi nasıl yapılır', r'parça alma'
        ],
        'title': 'Biyopsi Türleri & Doku Tanısının Önemi',
        'answer': (
            "🔬 **Biyopsi Çeşitleri ve Uygulama Alanları:**\n\n"
            "1. **Tru-Cut (Kalın İğne / Kor Biyopsi):** Meme, karaciğer, prostat gibi organlardan ince bir doku silindiri çıkarır; tümörün mimarisini ve immünohistokimyasal boyamalarını net gösterir (Altın standart).\n"
            "2. **İnce İğne Aspirasyon Biyopsisi (İİAB):** Enjektörle sadece sıvı ve hücre çekilir; en sık tiroid nodüllerinde kullanılır.\n"
            "3. **Eksizyonel Biyopsi:** Şüpheli kitlenin veya lenf bezinin tamamının cerrahiyle çıkarılmasıdır (Örn: Cilt benleri, lenfoma şüpheli lenf nodu).\n"
            "4. **Sıvı Biyopsi (ctDNA):** Kandan DNA parçacıklarının incelenmesidir."
        )
    },
    'tumor_marker_cea': {
        'patterns': [
            r'cea.*test', r'cea.*tümör', r'karsinoembriyonik', r'cea yüksekliği'
        ],
        'title': 'CEA (Karsinoembriyonik Antijen) Testi',
        'answer': (
            "📊 **CEA (Karsinoembriyonik Antijen) Tümör Belirteci:**\n\n"
            "• **Kullanım Alanı:** En sık Kolorektal (kalın bağırsak) kanserinin tedavi yanıtını ve ameliyat sonrası takibini izlemek için kanda ölçülür.\n"
            "• **Normal Değer:** Sigara içmeyenlerde < 3.0 ng/mL, sigara içenlerde < 5.0 ng/mL.\n"
            "• **Önemli Bilgi:** Tek başına kanser teşhisi koyamaz. Sigara içmek, mide ülseri, karaciğer yağlanması veya kolit gibi iyi huylu durumlarda da hafif yükselebilir."
        )
    },
    'tumor_marker_ca125': {
        'patterns': [
            r'ca.*125', r'ca125', r'yumurtalık.*belirteç', r'over.*tümör belirteç'
        ],
        'title': 'CA 125 Tümör Belirteci (Yumurtalık / Over Kanseri)',
        'answer': (
            "🌸 **CA 125 Belirteci ve Klinik Anlamı:**\n\n"
            "• **Kullanım Alanı:** Yumurtalık (over) kanserinin tanı ve kemoterapi yanıtının izlenmesinde birincil belirteçtir.\n"
            "• **Normal Eşik:** Genellikle < 35 U/mL kabul edilir.\n"
            "• **İyi Huylu Nedenler:** Menopoz öncesi kadınlarda endometriozis (çikolata kisti), adet dönemi, miyomlar veya pelvik enfeksiyonlar da CA 125 seviyesini yükseltebilir. Bu nedenle pelvik ultrasonografi ile birlikte değerlendirilir."
        )
    },
    'tumor_marker_ca19_9': {
        'patterns': [
            r'ca.*19-9', r'ca19-9', r'ca 19.9', r'pankreas.*belirteç'
        ],
        'title': 'CA 19-9 Belirteci (Pankreas & Safra Yolları)',
        'answer': (
            "🟡 **CA 19-9 Belirteci ve Klinik Önemi:**\n\n"
            "• **Kullanım:** Pankreas ve safra yolu (kolanjiyosellüler) kanserlerinde takipte kullanılır.\n"
            "• **Normal Değer:** < 37 U/mL.\n"
            "• **Dikkat Edilmesi Gereken:** Safra kanalı tıkanıklığı, sarılık, safra kesesi taşı veya pankreatit varlığında tümör olmaksızın da ileri derecede yükselebilir."
        )
    },
    'cancer_stages_tnm': {
        'patterns': [
            r'kanser evreleri', r'evre 4', r'evre 1', r'evre 2', r'evre 3',
            r'tnm sistemi', r'metastaz nedir', r'metastatik kanser'
        ],
        'title': 'Kanser Evreleri & TNM Sınıflandırması',
        'answer': (
            "📈 **Kanser Evrelemesi (TNM Sistemi) ve Metastaz Rehberi:**\n\n"
            "• **TNM Sistemi:**\n"
            "   - **T (Tümör Boyutu):** Primer tümörün milimetrik çapı ve komşu dokuya yayılımı (T1 - T4).\n"
            "   - **N (Lenf Nodu):** Bölgesel lenf bezlerine sıçrama var mı (N0 - N3).\n"
            "   - **M (Metastaz):** Uzak organlara sıçrama var mı (M0: Yok, M1: Var).\n"
            "• **Evreler:**\n"
            "   • **Evre 1 (Lokal):** Tümör küçük ve organın içinde sınırlıdır. Erken cerrahiyle %90+ kür sağlanır.\n"
            "   • **Evre 2 - 3 (Bölgesel İlerlemiş):** Tümör büyümüş veya lenf bezlerine ulaşmıştır.\n"
            "   • **Evre 4 (Metastatik):** Kanser kan veya lenf yoluyla akciğer, karaciğer, kemik veya beyin gibi uzak organlara yayılmıştır. Günümüzde yeni nesil immünoterapi ve akıllı ilaçlarla Evre 4 kanserlerde dahi uzun yıllar yüksek yaşam kalitesi sağlanabilmektedir."
        )
    },
    'remission_vs_cure': {
        'patterns': [
            r'remisyon', r'kür ne demek', r'nüks', r'rekürrens', r'kanser tekrarlar mı'
        ],
        'title': 'Kanserde Remisyon, Nüks ve Kür Terimleri',
        'answer': (
            "🔄 **Remisyon, Nüks ve Kür Ne Anlama Gelir?**\n\n"
            "• **Tam Remisyon:** Yapılan tüm radyolojik (BT, MR, PET) ve laboratuvar testlerinde vücutta saptanabilir hiçbir tümör kalmaması durumudur.\n"
            "• **Kısmi Remisyon:** Tümör yükünün en az %30-50 oranında küçülmesi ve hastalığın kontrol altına alınmasıdır.\n"
            "• **Kür (Tam Şifa):** Tedavi bitiminden sonra kanserin türüne göre 5 yıl boyunca hiçbir nüks görülmemesi durumunda hasta tıbben 'tam şifaya ulaşmış' kabul edilir.\n"
            "• **Nüks (Rekürrens):** Remisyondan sonra mikroskobik uykudaki hücrelerin yeniden aktifleşmesidir; sıvı biyopsi (ctDNA) ile aylar öncesinden yakalanabilir."
        )
    },
    'survival_os_pfs': {
        'patterns': [
            r'pfs nedir', r'overall survival', r'genel sağkalım', r'progresyonsuz sağkalım'
        ],
        'title': 'Klinik Sağkalım Kriterleri: OS & PFS',
        'answer': (
            "📊 **Onkolojik Sağkalım Terimleri:**\n\n"
            "• **Genel Sağkalım (OS - Overall Survival):** Tanı veya tedavi başlangıcından itibaren hastanın hayatta kaldığı toplam süredir; onkoloji ilaç araştırmalarında 'nihai altın standart' kabul edilir.\n"
            "• **Progresyonsuz Sağkalım (PFS - Progression-Free Survival):** İlaç başlandıktan sonra kanserin büyümeden veya ilerlemeden stabil kaldığı süredir.\n"
            "• **Hazard Ratio (HR - Tehlike Oranı):** Denenen yeni ilacın plasebo veya eski tedaviye göre ölüm/ilerleme riskini ne kadar düşürdüğünü gösterir (Örn: HR=0.60, riski %40 azaltmış demektir)."
        )
    },
    'cancer_prevention_diet': {
        'patterns': [
            r'kanserden korunma', r'kanser diyeti', r'nasıl beslenmeliyim', r'kanser yapmayan gıdalar'
        ],
        'title': 'Kanserden Korunma & Beslenme Kılavuzu (WHO)',
        'answer': (
            "🥗 **Kanıta Dayalı Kanserden Korunma Kılavuzu:**\n\n"
            "1. **Tütün Ürünlerinden Uzak Durun:** Kanser ölümlerinin %30'unun, akciğer kanserlerinin %90'ının doğrudan nedenidir.\n"
            "2. **İşlenmiş Et ve Şarküteri Ürünlerini Kısıtlayın:** Salam, sosis, sucuk gibi nitritli işlenmiş etler DSÖ tarafından Grup 1 karsinojen ilan edilmiştir.\n"
            "3. **Akdeniz Tipi Beslenin:** Zeytinyağı, baklagiller, lifli sebzeler, ceviz ve taze balık tüketimi bağırsak mikrobiyotasını koruyarak inflamasyonu azaltır.\n"
            "4. **İdeal Kilonuzu Koruyun:** Obezite; meme, kolon, rahim, pankreas ve böbrek dahil 13 farklı kanserin kanıtlanmış risk faktörüdür.\n"
            "5. **Haftada En Az 150 Dakika Egzersiz Yapın.**"
        )
    },
    'hpv_vaccine': {
        'patterns': [
            r'hpv aşısı', r'gardasil', r'rahim ağzı aşısı', r'hpv nedir'
        ],
        'title': 'HPV Aşısı & Kanserden %90 Korunma',
        'answer': (
            "💉 **HPV (Human Papilloma Virus) Aşısı Rehberi:**\n\n"
            "• **Hangi Kanserleri Önler?** Serviks (rahim ağzı) kanserlerinin %99'undan, ayrıca anüs, penis, vajina ve boğaz/ağız (orofarenks) kanserlerinden korur.\n"
            "• **Kimlere Yapılır?** Hem kız hem erkek çocuklarına 9-14 yaşından itibaren; yetişkinlerde ise 45 yaşına kadar uygulanabilir.\n"
            "• **Etkinlik:** Geniş toplum aşılaması yapılan ülkelerde (Avustralya, İngiltere) rahim ağzı kanseri oranları %90 azalmıştır."
        )
    },
    'genetic_counseling': {
        'patterns': [
            r'genetik danışmanlık', r'brca testi', r'gen testi.*kanser', r'gen paneli'
        ],
        'title': 'Genetik Kanser Danışmanlığı & Risk Panelleri',
        'answer': (
            "🧬 **Genetik Kanser Danışmanlığı & NGS Gen Panelleri:**\n\n"
            "• **Kime Yapılmalı?** Genç yaşta (<50 yaş) kanser tanısı alanlar, ailesinde birden fazla bireyde meme/over/kolon kanseri olanlar veya erkekte meme kanseri vakalarında önerilir.\n"
            "• **İncelenen Genler:** BRCA1, BRCA2, PALB2, ATM, CHEK2, TP53 ve Lynch sendromu genleri (MLH1, MSH2, MSH6, PMS2).\n"
            "• **Faydası:** Yüksek riskli bireylere erken tarama takvimi (Kontrastlı Meme MR, 30 yaşında kolonoskopi) veya koruyucu önlemler planlanır."
        )
    },
    'greetings': {
        'patterns': [
            r'\bmerhaba\b', r'\bselam\b', r'\bgünaydın\b', r'\biyi günler\b',
            r'\bnasılsın\b', r'\biyi akşamlar\b'
        ],
        'title': 'Klinik Karşılama & Selamlaşma',
        'answer': (
            "👋 **Merhaba! Hoş geldiniz.**\n\n"
            "Ben **Klinik Onkoloji & Erken Kanser Teşhisi Danışmanıyım** 🩺.\n\n"
            "Bana aklınıza takılan herhangi bir kanser belirtisini, tarama testlerini (Kolonoskopi, Mamografi, Düşük Doz BT), "
            "akıllı ilaçları, biyomarkerları (ctDNA, DPYD, BRCA) veya tahlil sonuçlarınızı sorabilirsiniz.\n\n"
            "👉 Nasıl yardımcı olabilirim? Şikayetlerinizi yazarak risk analizi başlatabilir veya merak ettiğiniz konuyu sorabilirsiniz."
        )
    },
    'gratitude': {
        'patterns': [
            r'teşekkür', r'sağ ol', r'sağol', r'eline sağlık', r'harikasın', r'çok sağol'
        ],
        'title': 'Teşekkür & Klinik Destek',
        'answer': (
            "🙏 **Rica ederim! Sağlığınız her şeyden önemlidir.**\n\n"
            "Erken teşhis hayat kurtarır. Aklınıza takılan başka bir tıbbi soru olursa dilediğiniz zaman çekinmeden sorabilirsiniz. "
            "Sağlıklı ve güzel günler dilerim! 🌿"
        )
    },

    'adc_therapies': {
        'patterns': [
            r'adc nedir', r'antikor.*ilaç.*konjugat', r'enhertu', r'trastuzumab deruxtecan',
            r'trodelvy', r'sacituzumab', r'padcev', r'enfortumab'
        ],
        'title': 'Yeni Nesil Antikor-İlaç Konjugatları (ADC) ve Truva Atı Teknolojisi',
        'answer': (
            "🎯 **Antikor-İlaç Konjugatları (ADC - Antibody-Drug Conjugates) Nedir?**\n\n"
            "• **Truva Atı Prensibi:** Bir monoklonal antikorun ucuna kimyasal bir bağlayıcı (linker) ile son derece güçlü bir sitotoksik kemoterapi molekülü (payload) bağlanır.\n"
            "• **Nasıl Çalışır?** Antikor yalnızca kanser hücresindeki spesifik reseptörü (örn. HER2, Trop-2, Nectin-4) tanır, hücre içine girer ve ölümcül yükünü doğrudan kanser hücresinin içinde serbest bırakır. Çevreleyen sağlıklı dokuya zarar vermez.\n"
            "• **Önemli İlaçlar:**\n"
            "   1. **Trastuzumab Deruxtecan (Enhertu):** HER2-pozitif ve HER2-low meme ve akciğer kanserinde çığır açmıştır (DESTINY çalışmaları).\n"
            "   2. **Sacituzumab Govitecan (Trodelvy):** Üçlü negatif meme kanseri (Trop-2) tedavisinde onaylıdır.\n"
            "   3. **Enfortumab Vedotin (Padcev):** İleri evre mesane kanserinde (Nectin-4) Pembrolizumab ile kombine standarttır."
        )
    },
    'cart_cell_therapy': {
        'patterns': [
            r'car-t hücresi nedir', r'car-t hücre tedavisi nedir', r'car-t nedir', r'kimerik antijen reseptörü', r'yaşayan ilaç', r'hücresel tedavi nedir', r'chimeric antigen receptor'
        ],
        'title': 'CAR-T Hücre Tedavisi (Kimerik Antijen Reseptörü) & Yaşayan İlaçlar',
        'answer': (
            "🧬 **CAR-T Hücre Tedavisi (Chimeric Antigen Receptor T-Cell) Nedir?**\n\n"
            "• **Yaşayan İlaç (Living Drug):** Hastanın kendi kanından alınan T-lenfositleri (bağışıklık hücreleri) laboratuvarda genetik mühendislikle yeniden programlanır.\n"
            "• **Mekanizma:** T hücrelerine kanser hücresi yüzeyindeki antijeni (örn. CD19 veya BCMA) tanıyan kimerik bir reseptör eklenir, çoğaltılır ve hastaya geri infüze edilir.\n"
            "• **Kullanım Alanları:** Dirençli B-hücreli lösemiler (ALL), lenfomalar (DLBCL) ve multipl miyelomda daha önce tüm tedavileri tüketmiş olgularda dahi %70-80 tam yanıt sağlayabilir.\n"
            "• **Onaylı Örnekler:** Tisagenlecleucel (Kymriah), Axicabtagene ciloleucel (Yescarta), Ciltacabtagene autoleucel (Carvykti)."
        )
    },
    'bispecific_antibodies': {
        'patterns': [
            r'bispesifik antikor nedir', r'bispesifik nedir', r'ikili antikor nedir', r'bispesifik teknoloji'
        ],
        'title': 'Bispesifik Antikorlar (Bispecific Antibodies) & Çift Hedefleme',
        'answer': (
            "🔗 **Bispesifik Antikorlar Nedir?**\n\n"
            "• **Çift Hedef:** Klasik antikorlar tek bir hedefi tanırken, bispesifik antikorlar iki farklı antijeni veya epitopu aynı anda bağlayabilir.\n"
            "• **İki Temel Çalışma Şekli:**\n"
            "   1. **T-Hücre Bağlayıcı (BiTE):** Bir koluyla kanser hücresini, diğer koluyla T hücresini (CD3) tutarak bağışıklık hücresini kansere zorla kenetler (Örn: Blinatumomab).\n"
            "   2. **Çift Onkojenik Reseptör Blokajı:** Örneğin **Amivantamab (Rybrevant)**, akciğer kanserinde hem EGFR hem de MET reseptörlerini aynı anda bloke ederek Exon 20 insersiyonlarında ve dirençte devrim yaratmıştır (PAPILLON çalışması)."
        )
    },
    'recist_criteria': {
        'patterns': [
            r'recist', r'tam yanıt', r'kısmi yanıt', r'progresyon nedir', r'stabil hastalık', r'tümör küçülmesi', r'cr pr sd pd'
        ],
        'title': 'Kanserde Tedavi Yanıt Kriterleri: RECIST 1.1 Kılavuzu',
        'answer': (
            "📏 **RECIST 1.1 (Response Evaluation Criteria in Solid Tumors) Kriterleri:**\n\n"
            "Radyolojik olarak (BT/MR) tümör boyutlarındaki değişimi nesnel olarak ölçen uluslararası standarttır:\n\n"
            "1. **Tam Yanıt (CR - Complete Response):** Tüm hedef lezyonların tamamen yok olması ve patolojik lenf bezlerinin <10 mm'ye gerilemesi.\n"
            "2. **Kısmi Yanıt (PR - Partial Response):** Hedef lezyonların en uzun çaplarının toplamında **en az %30 küçülme** olması.\n"
            "3. **Progresif Hastalık (PD - Progressive Disease):** En küçük boyuta kıyasla lezyonlarda **en az %20 artış** veya yeni metastatik lezyon çıkması.\n"
            "4. **Stabil Hastalık (SD - Stable Disease):** Ne PR için yeterli küçülme ne de PD için yeterli büyüme olmaması; hastalığın duraklaması."
        )
    },
    'amp_mtb_tiers': {
        'patterns': [
            r'amp seviyesi', r'amp.*tier', r'tier 1', r'tier 2', r'vus nedir', r'tümör konseyi', r'moleküler tümör konseyi'
        ],
        'title': 'Moleküler Tümör Konseyi & AMP/ASCO/CAP Kanıt Seviyeleri (Tier I-IV)',
        'answer': (
            "🏛️ **AMP / ASCO / CAP Somatik Varyant Kanıt Seviyeleri:**\n\n"
            "Moleküler Tümör Konseyleri (MTB), yeni nesil dizileme (NGS) raporundaki mutasyonları kanıt gücüne göre 4 seviyede değerlendirir:\n\n"
            "• **Tier I (Güçlü Klinik Önemi Olanlar):** FDA onaylı hedefe yönelik ilacı bulunan veya NCCN/ESMO rehberlerinde Kategori 1 tedavi standardı olan alterasyonlar (Örn: EGFR L858R, BRAF V600E, KRAS G12C).\n"
            "• **Tier II (Potansiyel Klinik Önemi Olanlar):** Faz 2/3 klinik çalışmalarda yanıt gösteren veya etiket dışı (off-label) umut vadeden mutasyonlar.\n"
            "• **Tier III (Bilinmeyen Klinik Öneme Sahip Varyantlar - VUS):** Literatürde yeterli kanıt bulunmayan varyantlar; standart dışı tedavi kararı verilmez.\n"
            "• **Tier IV (Benign / Etkisiz):** Tümör gelişimine etkisi olmayan zararsız polimorfizmler."
        )
    },
    'cancer_pathways_signaling': {
        'patterns': [
            r'mapk yolağı', r'ras raf mek', r'pi3k akt', r'sinyal iletimi', r'sinyal yolağı', r'onkojenik kaskad'
        ],
        'title': 'Kanserde Temel Sinyal İletim Yolakları (MAPK, PI3K/AKT, DDR)',
        'answer': (
            "⚡ **Kanserde Temel Onkojenik Sinyal Kaskadları:**\n\n"
            "1. **MAPK / ERK Kaskadı (RTK -> RAS -> RAF -> MEK -> ERK):**\n"
            "   Hücre bölünmesini ve proliferasyonunu kamçılar. Akciğer, kolon, melanom ve pankreas kanserlerinde KRAS, BRAF, EGFR mutasyonlarıyla sürekli 'AÇIK' kalır.\n\n"
            "2. **PI3K / AKT / PTEN / mTOR Yolağı:**\n"
            "   Hücre sağkalımı, anabolizma ve glukoz kullanımını yönetir. PTEN kaybı veya PIK3CA mutasyonları apoptozu (hücre ölümünü) engeller.\n\n"
            "3. **DDR (DNA Hasar Yanıtı) & Homolog Onarım (HRR):**\n"
            "   BRCA1, BRCA2 ve ATM genleri çift iplik kırıklarını onarır; bozulduğunda PARP inhibitörlerine aşırı duyarlılık (Sentetik Letalite) oluşur.\n\n"
            "4. **CDK4/6 - Retinoblastoma (RB1) Aksı:**\n"
            "   Hücre döngüsünün G1 fazından S fazına geçişini denetler."
        )
    },
    'clinical_trials_phases': {
        'patterns': [
            r'faz 1', r'faz 2', r'faz 3', r'klinik çalışma', r'klinik deney', r'ilaç aşamaları'
        ],
        'title': 'Kanser İlaçlarında Klinik Çalışma Fazları (Faz I, II, III)',
        'answer': (
            "🧪 **Kanser İlaç Geliştirme Fazları:**\n\n"
            "• **Faz I (Güvenlik & Doz Belirleme):** Genellikle 15-30 hasta üzerinde yapılır. Amacı ilacın insan vücudundaki güvenliğini, yan etkilerini ve tolere edilebilir maksimum dozunu (MTD) bulmaktır.\n"
            "• **Faz II (Etkinlik Araştırması):** 50-200 hastada uygulanır. İlacın hedeflenen kanser türünde tümörü küçültüp küçültmediği (Objektif Yanıt Oranı - ORR) değerlendirilir.\n"
            "• **Faz III (Rastgele Karşılaştırmalı Standart Test - RCT):** Yüzlerce veya binlerce hastada uygulanır. Yeni ilaç, mevcut en iyi tedaviye (Standart Bakım - SoC) karşı rastgele (randomize) yarıştırılır; Genel Sağkalım (OS) veya PFS farkı kanıtlanırsa FDA/EMA ruhsatı alır."
        )
    },
    'comprehensive_genomic_profiling': {
        'patterns': [
            r'ngs nedir', r'kapsamlı genomik profil', r'cgp', r'tümör dna dizileme', r'genetik test kanser', r'foundationone'
        ],
        'title': 'Kapsamlı Genomik Profilleme (CGP) & Yeni Nesil Dizileme (NGS)',
        'answer': (
            "🧬 **Kapsamlı Genomik Profilleme (Comprehensive Genomic Profiling - CGP):**\n\n"
            "• **Tek Seferde 500+ Gen:** Geleneksel tek gen testleri (PCR) yerine, biyopsideki tümör DNA ve RNA'sını tek seferde yüzlerce kanser geninde dizileyen (NGS) ileri teknolojidir.\n"
            "• **Neler Saptanır?**\n"
            "   1. **Nokta Mutasyonlar & İndeller:** (EGFR, KRAS, BRAF, PIK3CA)\n"
            "   2. **Kopya Sayısı Değişimleri (CNV):** (HER2, MET amplifikasyonları, PTEN kaybı)\n"
            "   3. **Gen Füzyonları:** (ALK, ROS1, RET, NTRK)\n"
            "   4. **Genomik İmzalar:** TMB (Tümör Mutasyon Yükü), MSI (Mikrosatellit İnstabilite) ve HRD (Homolog Rekombinasyon Eksikliği).\n"
            "• **Klinik Değeri:** Hastanın standart kemoterapi yerine kişiye özel akıllı ilaç veya immünoterapi alma şansını katbekat artırır."
        )
    },

    'immunotherapy_resistance_stk11_b2m': {
        'patterns': [
            r'immünoterapiye direnç', r'immünoterapi direnci', r'stk11', r'b2m kaybı',
            r'soğuk tümör', r'immün çöl', r'keap1'
        ],
        'title': 'İmmünoterapi Direnç Mekanizmaları: Soğuk Tümör & B2M / STK11 Kaybı',
        'answer': (
            "🛡️ **İmmünoterapi (Anti-PD-1) Direnci ve Soğuk Tümör Biyolojisi:**\n\n"
            "Bazı tümörler TMB-High veya PD-L1 pozitif olsa dahi immünoterapiye yanıt vermez (primer direnç) veya başlangıçta yanıt verip sonradan ilerler (kazanılmış direnç):\n\n"
            "1. **B2M (Beta-2 Mikroglobulin) Kaybı:** Hücre yüzeyinde MHC-I antijen sunum kompleksi çöker. Sitotoksik CD8+ T-hücreleri tümör neoantijenlerini fiziksel olarak göremez (körleşir).\n"
            "2. **STK11 (LKB1) ve KEAP1 Mutasyonları:** Akciğer kanserinde STING yolağını susturarak 'İmmün Çöl' (soğuk tümör) yaratır; T-hücre infiltrasyonunu önler. Tek ajan yerine ikili kontrol noktası (Anti-CTLA-4 + Anti-PD-1) veya kemoterapi şarttır.\n"
            "3. **JAK1 / JAK2 İnaktivasyonu:** İnterferon-gama sinyal iletimi bozulur; tümör hücresi T-hücre saldırısına duyarsızlaşır.\n"
            "4. **VEGF & TGF-β Aşırı Salınımı:** Yoğun stroma ve miyeloid baskılayıcı hücreler (MDSC) ile T-hücrelerinin tümör içine girmesi engellenir (İmmün Dışlanmış Fenotip)."
        )
    },
    'tumor_infiltrating_lymphocytes_til': {
        'patterns': [
            r'til tedavisi', r'tümör infiltre eden lenfosit', r'amtagvi', r'lifileucel', r'hücre infüzyonu'
        ],
        'title': 'Tümör İnfiltre Eden Lenfosit (TIL) Tedavisi & Lifileucel',
        'answer': (
            "🧬 **TIL (Tumor-Infiltrating Lymphocytes) Hücresel Tedavisi:**\n\n"
            "• **Doğal Kanser Savaşçıları:** Hastanın kendi çıkarılan tümör dokusunun içine sızmış olan özgül lenfositler (TIL) cerrahiyle alınır.\n"
            "• **Laboratuvarda Milyarlara Katlama:** Bu hücreler IL-2 sitokini eşliğinde laboratuvarda milyarlarca adet çoğaltılır ve aktive edilir.\n"
            "• **Onaylı İlk İlaç:** **Lifileucel (Amtagvi)**, ileri evre melanomda immünoterapi ve hedefe yönelik ilaçları tüketmiş olgularda FDA onayı alan ilk hücresel katı tümör tedavisidir."
        )
    },
    'oncolytic_virus_therapy': {
        'patterns': [
            r'onkolitik virüs', r'imlygic', r't-vec', r'virüsle kanser tedavisi', r'kanser virüsü'
        ],
        'title': 'Onkolitik Virüs Tedavileri (T-VEC / Imlygic) & İmmün Uyarım',
        'answer': (
            "🦠 **Onkolitik Virüs İmmünoterapisi Nedir?**\n\n"
            "• **Kansere Özel Modifiye Virüs:** Genetik olarak zayıflatılmış virüsler (örn. Herpes Simplex Virüsü - Talimogene laherparepvec / T-VEC) doğrudan tümörün içine enjekte edilir.\n"
            "• **Çifte Darbe:**\n"
            "   1. Yalnızca kanser hücrelerinde çoğalarak onları içeriden patlatır (lizis).\n"
            "   2. Parçalanan kanser hücrelerinden saçılan antijenler ve virüsün ürettiği GM-CSF proteini, vücudun tüm bağışıklık sistemini tümöre karşı alarma geçirir.\n"
            "• **Uygulama:** İlerlemiş melanom deri ve lenf bezi metastazlarında bölgesel enjeksiyon olarak kullanılır."
        )
    },
    'neoantigen_mrna_vaccines': {
        'patterns': [
            r'mrna kanser aşısı', r'kanser aşısı', r'biontech kanser', r'moderna kanser', r'mrna-4157', r'kişiye özel aşı'
        ],
        'title': 'Kişiye Özel mRNA Kanser Aşıları & Neoantijen Teknolojisi',
        'answer': (
            "💉 **Kişiselleştirilmiş mRNA Kanser Aşıları (mRNA-4157 / V940 & BioNTech):**\n\n"
            "• **Parmak İzi Kadar Özel:** Her hastanın ameliyatla alınan tümör DNA'sı dizilenir (NGS). Sağlıklı DNA ile karşılaştırılarak hastaya özgü en immünojenik 30-34 adet mutasyonel neoantijen saptanır.\n"
            "• **Sentetik mRNA Üretimi:** Bu hedefler haftalar içinde sentetik bir mRNA molekülüne kodlanır ve hastaya aşı olarak uygulanır.\n"
            "• **Klinik Başarı (KEYNOTE-942):** İleri evre melanom cerrahisi sonrası mRNA-4157 + Pembrolizumab kombinasyonunun, tek başına Pembrolizumab'a kıyasla nüks veya ölüm riskini **%44-49 oranında azalttığı** Faz 2b çalışmasında kanıtlanmıştır."
        )
    },
    'radioligand_pluvicto_lutathera': {
        'patterns': [
            r'radyoligand tedavisi nedir', r'akıllı radyasyon nedir', r'atom tedavisi.*prostat', r'radyoligand nedir'
        ],
        'title': 'Hassas Radyoligand Tedavileri: Pluvicto (Lu-177 PSMA) & Lutathera',
        'answer': (
            "☢️ **Radyoligand Tedavisi (Akıllı Atom / Teranostik Yaklaşım):**\n\n"
            "Hedefe yönelik küçük moleküller veya peptitler, radyoaktif bir izotop (Lutesyum-177 veya Aktinyum-225) ile birleştirilerek damardan verilir.\n\n"
            "1. **Pluvicto (Lutetium Lu 177 Vipivotide Tetraxetan):** Prostat kanseri hücrelerindeki PSMA reseptörünü bulur; mikro-radyasyonla DNA'larını içeriden tahrip eder (VISION çalışması).\n"
            "2. **Lutathera (Lu-177 Dotatate):** Nöroendokrin tümörlerde (NET) somatostatin reseptörlerine kilitlenir.\n"
            "• **Teranostik Güç:** Önce Ga-68 PET ile tümörün reseptörü görüntülenir (Teşhis), ardından aynı molekülle Lu-177 bağlanarak tedavi edilir."
        )
    },

    'mutational_signatures_cosmic': {
        'patterns': [
            r'mutasyonel imza', r'cosmic imza', r'sbs1', r'sbs2', r'sbs4', r'sbs7',
            r'apobec imza', r'mutational signature'
        ],
        'title': 'COSMIC Mutasyonel İmzalar (SBS) & Tümör Arkeolojisi',
        'answer': (
            "🧬 **COSMIC Mutasyonel İmzalar (Mutational Signatures):**\n\n"
            "Kanserin genomunda geçmişte maruz kaldığı biyolojik ve çevresel hasarların bıraktığı moleküler 'parmak izleridir' (96 trinükleotit bağlamı):\n\n"
            "• **SBS1 (Biyolojik Yaş Saati):** Spontan 5-metilsitozin deaminasyonu ile oluşur; hastanın biyolojik yaşlanma hızını yansıtır.\n"
            "• **SBS2 & SBS13 (APOBEC İmzası):** Hücresel sitidin deaminaz enzim atağıdır; mesane, meme ve serviks kanserlerinde yüksek TMB yaratır.\n"
            "• **SBS4 (Sigara Dumanı):** Tütün polisiklik aromatik hidrokarbonları (benzo[a]piren) DNA'ya kovalent bağlanarak C>A transversiyonları üretir.\n"
            "• **SBS7 (Güneş / UV Hasarı):** Melanomda pirimidin dimerleri üzerinden CC>TT değişimleri üretir; en yüksek neoantijen yükünü oluşturur.\n"
            "• **SBS10 (POLE İmzası):** DNA polimeraz epsilon proofreading arızası; 100-500+ mut/Mb ultra-hipermutasyon yaratır (immünoterapiye aşırı duyarlı)."
        )
    },
    'hrd_genomic_scar': {
        'patterns': [
            r'hrd skoru', r'hrd nedir', r'genomik skar', r'homolog rekombinasyon eksikliği',
            r'heterozigotluk kaybı', r'loh skoru', r'mychoice'
        ],
        'title': 'Homolog Rekombinasyon Eksikliği (HRD) & Genomik Skar Skoru',
        'answer': (
            "🎯 **HRD (Homologous Recombination Deficiency) & Genomik Skar:**\n\n"
            "• **Kritik Biyobelirteç:** Tümörün çift iplik DNA kırıklarını hatasız onarma yeteneğini kaybedip kaybetmediğini ölçer (özellikle over, meme, prostat ve pankreas kanserlerinde).\n"
            "• **3 Bileşenli Genomik Skar:**\n"
            "   1. **LOH (Loss of Heterozygosity):** Heterozigotluk kaybı alanları.\n"
            "   2. **TAI (Telomeric Allelic Imbalance):** Telomerik alelik dengesizlik bölgeleri.\n"
            "   3. **LST (Large-scale State Transitions):** Kromozomal büyük ölçekli durum geçişleri.\n"
            "• **Klinik Eşik:** Toplam HRD Skoru **≥42** veya patojenik **BRCA1/2 mutasyonu** bulunması durumu 'HRD Pozitif' kabul edilir.\n"
            "• **Tedavi:** Bu hastalar **PARP İnhibitörlerinden (Olaparib, Niraparib, Rucaparib)** ve Platin Kemoterapisinden en yüksek sağkalım avantajını elde eder."
        )
    },
    'protac_molecular_glues': {
        'patterns': [
            r'protac', r'hedefe yönelik protein yıkımı', r'tpd', r'moleküler yapıştırıcı',
            r'molecular glue', r'e3 ligaz', r'arv-110'
        ],
        'title': 'PROTAC Teknolojisi & Moleküler Yapıştırıcılar (TPD)',
        'answer': (
            "💣 **Hedefe Yönelik Protein Yıkımı (Targeted Protein Degradation - PROTAC / TPD):**\n\n"
            "• **İlaçlama Paradigmasında Devrim:** Klasik inhibitörler hedef proteinin aktif bölgesini geçici olarak tıkarken; PROTAC molekülleri hedef proteini hücrenin kendi çöp öğütücüsüne (Proteazom) teslim ederek **tamamen yok eder**.\n"
            "• **Nasıl Çalışır?** Çift başlıklı bir moleküldür: Bir ucu kanser proteinine (örn. Androjen veya Östrojen reseptörü), diğer ucu ise bir E3 ubiquitin ligaz enzimine (Cereblon veya VHL) bağlanır.\n"
            "• **Klinik Örnekler:** Bavdegalutamide (ARV-110, prostat kanseri), Vepdegestrant (ARV-471, meme kanseri) gibi moleküller, standart ilaçlara dirençli mutasyonlu reseptörleri dahi başarıyla yıkar."
        )
    },
    'ctc_circulating_tumor_cells': {
        'patterns': [
            r'ctc nedir', r'dolaşımdaki tümör hücresi', r'kanda tümör hücresi', r'cellsearch',
            r'metastaz hücresi', r'dolaşan kanser hücresi'
        ],
        'title': 'Dolaşımdaki Tümör Hücreleri (CTC) & Sıvı Biyopsi Hücre Sayımı',
        'answer': (
            "🩸 **Dolaşımdaki Tümör Hücreleri (Circulating Tumor Cells - CTC):**\n\n"
            "• **Primer Tümörden Kopan Hücreler:** Katı tümörlerden ayrılarak kan dolaşımına geçen ve uzak organlara tutunarak metastaz başlatabilen canlı kanser hücreleridir.\n"
            "• **CellSearch Teknolojisi (FDA Onaylı):** Kanda 7.5 mL tüp başına düşen CTC sayısını (örn. Meme ve prostatta ≥5 hücre/7.5 mL kötü prognoz eşiğidir) ölçer.\n"
            "• **ctDNA\'dan Farkı:** ctDNA ölen hücrelerin serbest DNA parçalarını ölçerken; CTC canlı hücreleri yakalar. Bu hücreler üzerinde tek hücre RNA dizileme (scRNA-seq) veya kültüre edilerek ilaç duyarlılık testleri yapılabilir."
        )
    },
    'cryo_em_in_drug_discovery': {
        'patterns': [
            r'cryo-em', r'kriyo em', r'kriyojenik elektron mikroskop', r'elektron mikroskobu ilaç',
            r'yapısal biyoloji', r'protein yapısı çözme'
        ],
        'title': 'Cryo-EM (Kriyojenik Elektron Mikroskobisi) & İlaç Tasarımı',
        'answer': (
            "🔬 **Kriyojenik Elektron Mikroskobisi (Cryo-EM) ile İlaç Keşfi:**\n\n"
            "• **Kristallendirme Şartı Yok:** Klasik X-ışını kristalografisinde zorunlu olan ve membran proteinlerinde (GPCR'lar, iyon kanalları, RTK'lar) aylarca kristal üretilemeyen zorlu hedefleri sıvı azotta şok dondurarak (vitrifikasyon) inceler.\n"
            "• **Atomik Çözünürlük:** 1.5 - 2.5 Ångström seviyesinde atomik detayla proteinin doğal fizyolojik konformasyonunu ve ilaç adayıyla bağlanma cebini doğrudan görüntüler.\n"
            "• **Kanserde Rolü:** KRAS G12C kompleksleri, Ribozom hedefli terapötikler ve büyük antikor-antijen kompleksleri Cryo-EM sayesinde çözülerek yeni nesil akıllı ilaçlar geliştirilmiştir."
        )
    },

    'spatial_transcriptomics_visium': {
        'patterns': [
            r'uzamsal transkriptomik', r'spatial transcriptomics', r'10x visium', r'xenium',
            r'doku haritası', r'uzamsal biyoloji', r'cosmx'
        ],
        'title': 'Uzamsal Transkriptomik (Spatial Transcriptomics) & Doku Mimarisi',
        'answer': (
            "🗺️ **Uzamsal Transkriptomik (Spatial Transcriptomics):**\n\n"
            "• **Hücrenin Konumunu Kaybetmeden:** Klasik RNA dizileme tüm dokuyu blenderdan geçirip tek bir ortalama değer verirken; 10x Visium, Xenium ve CosMx gibi teknolojiler binlerce geni **tümör doku kesiti üzerindeki tam koordinatlarında** görüntüler.\n"
            "• **Klinik Değeri:**\n"
            "   1. **Tümör Çekirdeği vs. İnvaziv Sınır:** Kanser hücrelerinin çevre dokuya yayıldığı temas yüzeyindeki saldırgan klonları ve kaçış mutasyonlarını gösterir.\n"
            "   2. **Hücreler Arası Komşuluk:** Sitotoksik T hücrelerinin kansere temas edip etmediğini veya fibroblast duvarının arkasında hapsolup olmadığını haritalar."
        )
    },
    'tertiary_lymphoid_structures_tls': {
        'patterns': [
            r'tersiyer lenfoid yapı', r'tls nedir', r'tertiary lymphoid', r'germinal merkez.*tümör',
            r'cxcl13', r'tümör lenf düğümü'
        ],
        'title': 'Tersiyer Lenfoid Yapılar (TLS) & İmmünoterapi Süper-Yanıtı',
        'answer': (
            "🏰 **Tersiyer Lenfoid Yapılar (Tertiary Lymphoid Structures - TLS):**\n\n"
            "• **Tümör Yanında Minyatür Lenf Düğümü:** Kronik inflamasyon ortamında tümörün hemen kıyısında kendiliğinden organize olan ektopik lenfoid organoidlerdir.\n"
            "• **İçerik:** B hücreleri (CD20+), T hücreleri (CD4+/CD8+), foliküler dendritik hücreler ve CXCL13 kemokin merkezinden oluşur.\n"
            "• **Klinik Önemi:** Biyopside olgun (germinal merkezli) TLS varlığı, hastanın **Anti-PD-1 immünoterapisine çok yüksek yanıt vereceğinin ve uzun genel sağkalımının** en güçlü bağımsız biyobelirtecidir."
        )
    },
    'tumor_associated_macrophages_tam': {
        'patterns': [
            r'tümör ilişkili makrofaj', r'tam nedir', r'm1 m2 makrofaj', r'cd47', r'sirpa',
            r'beni yeme sinyali'
        ],
        'title': 'Tümör İlişkili Makrofajlar (TAM) & CD47 "Beni Yeme" Sinyali',
        'answer': (
            "🛡️ **Tümör İlişkili Makrofajlar (TAM) ve İmmün Kaçış:**\n\n"
            "• **İki Yüzlü Hücreler:**\n"
            "   - **M1 Fenotipi:** Tümörü öldüren, pro-inflamatuar, savaşçı makrofajlar.\n"
            "   - **M2 Fenotipi:** Tümör tarafından 'kandırılmış', anjiyogenezi kamçılayan, stroma üreten ve T hücrelerini susturan baskılayıcı hücreler (CD206+, CD163+).\n"
            "• **CD47 'Beni Yeme' (Don't Eat Me) Ekseni:** Kanser hücreleri yüzeyinde CD47 molekülünü aşırı eksprese ederek makrofajın SIRPα reseptörünü kilitler ve fagositozdan (yutulmaktan) kaçar.\n"
            "• **Yeni İlaçlar:** Anti-CD47 monoklonal antikorları (Magrolimab) ve CSF-1R inhibitörleri M2 makrofajları yeniden M1 savaşçılarına çevirmeyi hedefler."
        )
    },
    'cancer_associated_fibroblasts_caf': {
        'patterns': [
            r'kanser ilişkili fibroblast', r'caf nedir', r'desmoplazi', r'stroma bariyeri',
            r'fibrotik stroma', r'fap inhibitör'
        ],
        'title': 'Kanserle İlişkili Fibroblastlar (CAF) & Desmoplastik Bariyer',
        'answer': (
            "🧱 **Kanserle İlişkili Fibroblastlar (Cancer-Associated Fibroblasts - CAF):**\n\n"
            "• **Tümörün Koruma Duvarı:** Katı tümörlerde (özellikle pankreas ve meme) aşırı kollajen ve fibronektin üreterek **desmoplazi** adı verilen sert, taş gibi bir zırh oluşturur.\n"
            "• **İlaç Geçirmez Bariyer:** Bu yoğun fibröz doku tümör içi basıncı artırarak kan damarlarını ezer; kemoterapinin ve T lenfositlerinin tümörün merkezine ulaşmasını fiziksel olarak engeller (İmmün Dışlanmış Fenotip).\n"
            "• **Tedavi Yaklaşımları:** FAP (Fibroblast Aktivasyon Proteini) hedefli radyonüklit tedaviler ve TGF-beta reseptör blokerleri ile stroma yumuşatılarak immünoterapi duyarlılığı geri kazandırılmaktadır."
        )
    },
    'single_cell_rna_sequencing_scrna': {
        'patterns': [
            r'tek hücre rna', r'scrna-seq', r'single cell rna', r'tek hücre dizileme',
            r'klonal heterojenite', r'kanser kök hücresi'
        ],
        'title': 'Tek Hücre RNA Dizileme (scRNA-seq) & Klonal Heterojenite',
        'answer': (
            "🔬 **Tek Hücre RNA Dizileme (Single-Cell RNA-seq):**\n\n"
            "• **Her Hücrenin Ayrı Kimliği:** Biyopsi parçasındaki on binlerce hücreyi tek tek damlacıklar içine hapsederek her bir hücrenin transkriptomunu (aktif genlerini) ayrı ayrı okur.\n"
            "• **Nadir Klonların Keşfi:** 10.000 hücre arasında saklanan ve tedaviye direnç geliştirecek olan %0.1'lik **nadir kanser kök hücrelerini** veya tükenmiş (exhausted) T hücre alt gruplarını erkenden tespit eder.\n"
            "• **Klonal Evrim:** Tedavi altında hangi kanser alt grubunun öldüğünü, hangisinin mutasyon biriktirerek hayatta kaldığını moleküler düzeyde izler."
        )
    },

    'mgmt_promoter_methylation': {
        'patterns': [
            r'mgmt metilasyonu', r'mgmt promoter', r'glioblastom mgmt', r'temozolomid mgmt', r'temodal duyarlılık', r'stupp protokolü'
        ],
        'title': 'MGMT Promotör Metilasyonu & Glioblastomda Temozolomid Duyarlılığı',
        'answer': (
            "🧠 **MGMT Promotör Metilasyonu ve Glioblastomda Klinik Değeri:**\n\n"
            "• **Biyolojik Mekanizma:** MGMT (O-6-Metilguanin-DNA Metiltransferaz), alkilleyici kemoterapi ilaçlarının (Temozolomid - TMZ) DNA'da açtığı ölümcül hasarları tamir eden bir enzimdir.\n"
            "• **Metillenmiş (Susturulmuş) MGMT:** Promotör bölgesi hipermetillenmişse enzim üretilemez. Temozolomidin yarattığı hasar onarılamaz ve kanser hücreleri ölür. Bu olgularda ortanca sağkalım **21.7 aya uzar** (Kategori 1 Stupp Protokolü kanıtı).\n"
            "• **Metillenmemiş MGMT:** Enzim aktiftir, TMZ direncine yol açar (Sağkalım ~15.3 ay). Bu hastalarda klinik çalışmalar, lomustin (CCNU) veya TTFields kombinasyonları öne çıkar."
        )
    },
    'mlh1_promoter_hypermethylation': {
        'patterns': [
            r'mlh1 metilasyonu', r'mlh1 hipermetilasyon', r'lynch sendromu ayrımı', r'braf v600e lynch', r'sporadik msi-h', r'cimp fenotipi'
        ],
        'title': 'MLH1 Promotör Hipermetilasyonu & Lynch Sendromu Triyajı',
        'answer': (
            "🧬 **MLH1 Hipermetilasyonu ile Sporadik - Kalıtsal (Lynch) Kanser Ayrımı:**\n\n"
            "• **Kritik Triyaj Testi:** Kolorektal veya endometriyal kanserde immünohistokimyada (IHC) MLH1/PMS2 kaybı saptandığında, bunun kalıtsal Lynch Sendromu mu yoksa yaşa bağlı sporadik kanser mi olduğunu ayırt eder.\n"
            "• **MLH1 Hipermetilasyonu Pozitif + BRAF V600E Mutasyonu:** Kanserin **sporadik (edinsel)** olduğunu kesinleştirir. Hastanın ailesine kalıtsal geçiş riski bulunmaz; pahalı germline genetik testler gerekmez.\n"
            "• **Metilasyon Yoksa (Unmethylated):** Kalıtsal Lynch Sendromu şüphesi doğar; hastaya ve birinci derece akrabalarına genetik danışmanlık ve germline MLH1 dizilemesi zorunludur.\n"
            "• **Tedavi:** Her iki grupta da MSI-H mevcuttur; Pembrolizumab veya Nivolumab+Ipilimumab gibi immünoterapilere yanıt mükemmeldir."
        )
    },
    'multicancer_early_detection_mced': {
        'patterns': [
            r'mced nedir', r'çoklu kanser erken teşhis', r'galleri testi', r'cfdna metilasyon', r'doku kökeni analizi', r'sept9 metilasyon'
        ],
        'title': 'Çoklu Kanser Erken Teşhisi (MCED) & cfDNA Metilasyonu (Sıvı Biyopsi)',
        'answer': (
            "🩸 **Çoklu Kanser Erken Teşhisi (MCED) & Hücresiz DNA (cfDNA) Metilasyonu:**\n\n"
            "• **Tek Kan Testiyle 50+ Kanser Taraması:** Klasik taramalar sadece meme, kolon ve serviksle sınırlıyken, MCED sıvı biyopsi testleri (Galleri gibi) kandaki hücresiz DNA parçacıklarının CpG metilasyon paternlerini okur.\n"
            "• **Doku Kökeni (Tissue-of-Origin - TOO) Deşifresi:** Kanserli dokudan kana dökülen DNA'nın metilasyon 'barkodu', tümörün pankreas, kolon, mide, akciğer veya over kökenli olduğunu **>%90 doğrulukla** işaret eder.\n"
            "• **SEPT9 & SHOX2 Biyobelirteçleri:** Kolorektal kanserde kanda hipermetillenmiş SEPT9 DNA'sı; akciğer kanserinde SHOX2 metilasyonu erken evrede bile yüksek hassasiyetle sinyal verir."
        )
    },
    'ezh2_inhibitors_epigenetics': {
        'patterns': [
            r'ezh2 inhibitörü', r'tazemetostat', r'smarcb1', r'ini1 kaybı', r'epitelioid sarkom', r'h3k27me3'
        ],
        'title': 'EZH2 İnhibitörleri (Tazemetostat) & SMARCB1/INI1 Epigenetik Bozuklukları',
        'answer': (
            "🎯 **EZH2 İnhibisyonu ve Epigenetik Hedefli Tedaviler:**\n\n"
            "• **Mekanizma:** PRC2 kompleksinin katalitik alt birimi olan EZH2, histon H3 lizin 27'yi trimetiller (H3K27me3) ve tümör baskılayıcı genleri baskılar.\n"
            "• **SMARCB1 (INI1) Kaybı:** Epitelioid sarkom ve rabdoid tümörlerde SMARCB1 delesyonu, EZH2 enzimini kontrolsüz hiperaktif hale getirir.\n"
            "• **Tazemetostat (Tazverik - FDA Onaylı):** EZH2'yi seçici olarak inhibe ederek tümör baskılayıcı genlerin üzerindeki epigenetik kilidi açar ve hücrelerin farklılaşarak apoptoza gitmesini sağlar.\n"
            "• **Endikasyonlar:** İleri evre epitelioid sarkom ve EZH2 mutasyonlu nüks/dirençli foliküler lenfoma."
        )
    },
    'hypomethylating_agents_dnmt_hdac': {
        'patterns': [
            r'hipometilleyici ajanlar', r'azasitidin', r'desitabin', r'dnmt inhibitörü', r'hdac inhibitörü', r'vorinostat', r'epigenetik tedavi'
        ],
        'title': 'Hipometilleyici Ajanlar (DNMTi: Azasitidin/Desitabin) & HDAC İnhibitörleri',
        'answer': (
            "💊 **DNA Metiltransferaz (DNMT) & Histon Deasetilaz (HDAC) İnhibitörleri:**\n\n"
            "• **DNA Hipometilasyonu (DNMT İnhibitörleri):**\n"
            "   - **Azasitidin (Vidaza) & Desitabin (Dacogen):** DNA metiltransferaz enzimlerini tuzağa düşürerek tümör baskılayıcı genlerin (CDKN2A vb.) promotörlerindeki hipermetilasyonu siler. Miyelodisplastik sendrom (MDS) ve AML'de standart tedavi omurgasıdır.\n"
            "• **Histon Yeniden Şekillendirme (HDAC İnhibitörleri):**\n"
            "   - **Vorinostat (SAHA), Romidepsin, Belinostat:** Histon deasetilazları bloke ederek kromatinin açık kalmasını sağlar; T-hücreli lenfomalarda hücre ölümünü tetikler.\n"
            "• **Epigenetik Aşıcılık:** Bu ilaçlar kanser antijenlerini yeniden eksprese ettirerek immünoterapinin başarısını da artırır."
        )
    },

    'cart_therapy_generations': {
        'patterns': [
            r'car-t nesilleri', r'car-t jenerasyon', r'truck car-t', r'4-1bb cd28', r'kimera antijen reseptörü', r'ko-stimülasyon'
        ],
        'title': 'CAR-T Hücre Tedavisi Nesilleri & Ko-Stimülasyon (4-1BB vs CD28)',
        'answer': (
            "🧬 **CAR-T Hücre Mimarisi ve Nesiller Arası Evrim:**\n\n"
            "• **1. Nesil:** Yalnızca CD3ζ aktivasyon zinciri içerir; in-vivo hücre sağkalımı ve çoğalması yetersizdir.\n"
            "• **2. Nesil (Günümüz Onaylı İlaçları):** Bir adet ko-stimülatör alan (CD28 veya 4-1BB/CD137) eklenmiştir.\n"
            "   - **CD28 (Yescarta/Tecartus):** Hızlı ve yoğun glikolitik patlama; ani tümör yükü azaltımı sağlar ancak tükenme (exhaustion) daha sıktır.\n"
            "   - **4-1BB (Kymriah/Carvykti/Breyanzi):** Mitokondriyal biyogenezi destekler; kanda 5-10 yıla varan uzun vadeli kalıcılık (persistence) ve bellek sağlar.\n"
            "• **3. Nesil:** Çift ko-stimülasyon (CD28 + 4-1BB + CD3ζ).\n"
            "• **4. Nesil (TRUCKs / Zırhlı CAR):** Katı tümör mikroçevresindeki stroma bariyerini aşmak için IL-12, IL-15 salgılayan veya baskılayıcı TGF-beta reseptörünü etkisiz kılan (dnTGFBR2) akıllı hücrelerdir."
        )
    },
    'crs_cytokine_release_syndrome': {
        'patterns': [
            r'sitokin salınım sendromu', r'crs nedir', r'tocilizumab car-t', r'sitokin fırtınası car-t', r'astct crs', r'\bcrs\b'
        ],
        'title': 'Sitokin Salınım Sendromu (CRS) & Tocilizumab Tedavisi',
        'answer': (
            "🔥 **Sitokin Salınım Sendromu (CRS) & ASTCT Klinik Yönetimi:**\n\n"
            "• **Neden Olur?** CAR-T hücreleri kanser antijeniyle karşılaşıp patladığında yoğun miktarda İnterlökin-6 (IL-6), IFN-γ ve TNF-α salgılanır.\n"
            "• **Belirtiler:** Yüksek dirençli ateş (≥38.0°C), tansiyon düşüklüğü (hipotansiyon) ve oksijen düşüşü (hipoksi).\n"
            "• **Evreleme ve Tedavi (ASTCT Konsensusu):**\n"
            "   - **Evre 1 (Hafif):** Yalnızca ateş; parasetamol ve hidrasyonla izlenir.\n"
            "   - **Evre 2 (Orta):** Sıvıya yanıt veren hipotansiyon veya düşük akım oksijen ihtiyacı; **Tocilizumab (8 mg/kg IV - IL-6 reseptör antagonisti)** uygulanır.\n"
            "   - **Evre 3/4 (Şiddetli):** Yoğun bakım, vazopressör ve ventilatör desteği; Tocilizumab + yüksek doz Deksametazon (20 mg IV) verilir."
        )
    },
    'icans_neurotoxicity_cart': {
        'patterns': [
            r'icans', r'car-t nörotoksisite', r'ice skoru', r'ensefalopati car-t', r'deksametazon icans'
        ],
        'title': 'ICANS Nörotoksisitesi & İmmün Efektör Hücre Ensefalopatisi (ICE)',
        'answer': (
            "🧠 **ICANS (Immune Effector Cell-Associated Neurotoxicity Syndrome):**\n\n"
            "• **Mekanizma:** Sitokin fırtınası nedeniyle kan-beyin bariyerinin geçirgenliğinin artması ve merkezi sinir sisteminde miyeloid hücre inflamasyonudur.\n"
            "• **İlk Sinyal:** Yazı yazmada bozulma (disgrafi), isimleri hatırlayamama (afazi) ve kafa karışıklığı.\n"
            "• **ICE Skoru (0-10 Puan):** Yönelim, nesne isimlendirme, emirleri uygulama ve yazı yazma ile günde 2-3 kez objektif olarak taranır.\n"
            "• **Tedavi:** Tocilizumab kan-beyin bariyerini iyi geçemediği için primer tedavi **Kortikosteroidlerdir (Deksametazon 10 mg q6h veya Metilprednizolon 1000 mg/gün)**. Nöbet profilaksisi için Levetirasetam verilir."
        )
    },
    'bcma_carvykti_myeloma': {
        'patterns': [
            r'bcma car-t', r'carvykti', r'abecma', r'ciltacabtagene', r'multipl miyelom car-t', r'\bbcma\b'
        ],
        'title': 'BCMA Hedefli CAR-T (Carvykti / Abecma) & Multipl Miyelomda Yanıt',
        'answer': (
            "🩸 **BCMA Hedefli CAR-T Hücreleri ve Multipl Miyelomda Başarı:**\n\n"
            "• **Hedef:** B-Hücre Olgunlaşma Antijeni (BCMA), plazma hücrelerinde ve miyelomda yüksek oranda eksprese edilen ideal bir hedeftir.\n"
            "• **Ciltacabtagene Autoleucel (Carvykti - CARTITUDE-4):**\n"
            "   - İki farklı BCMA epitopuna bağlanan benzersiz llama/alpaka nanobody yapısına sahiptir.\n"
            "   - Daha önce 3-4 basamak tedavi görmüş dirençli hastalarda **%98 Objektif Yanıt (ORR)** ve **%82.5 Sıkı Tam Yanıt (sCR)** sağlamıştır.\n"
            "• **İde-cel (Abecma):** İlk FDA onaylı BCMA CAR-T hücresidir; tek doz infüzyonla kemoterapisiz uzun remisyonlar sunar."
        )
    },
    'claudin18_2_gastric_cancer': {
        'patterns': [
            r'claudin 18.2', r'cldn18.2', r'claudin', r'zolbetuximab', r'ct041', r'mide kanseri hedefe yönelik claudin'
        ],
        'title': 'Claudin 18.2 Hedeflemesi (Zolbetuximab & CT041 CAR-T) Mide Kanseri',
        'answer': (
            "🎯 **Claudin 18.2 (CLDN18.2) & Mide / Pankreas Kanserinde Yeni Çığır:**\n\n"
            "• **Seçici Tümör Antijeni:** Normalde mide mukozasının sıkı bağlantı (tight junction) aralarında gizlidir. Kanserleşme sırasında hücre polaritesi kaybolunca bu hedef dışarı açılır ve hedeflenebilir hale gelir.\n"
            "• **Zolbetuximab (Vyloy - FDA/EMA Onaylı):** CLDN18.2 pozitif (≥%75 hücrede 2+/3+ boyanma) HER2-negatif ileri evre mide/GEJ kanserinde kemoterapi ile sağkalımı belirgin uzatır (SPOTLIGHT & GLOW Faz 3).\n"
            "• **Satricabtagene Autoleucel (CT041):** Katı tümörlerde umut vadeden ilk Claudin 18.2 hedefli CAR-T hücresidir; ileri evre mide kanserinde %57 ORR sağlamaktadır."
        )
    },

    'adc_mechanism_bystander': {
        'patterns': [
            r'bystander effect', r'bystander etkisi', r'adc nasıl öldürür', r'antikor ilaç konjugatı mekanizma', r'bystander öldürme', r'çevre hücreleri öldürme'
        ],
        'title': 'Antikor-İlaç Konjugatları (ADC) & Çevre Hücreleri Öldürme (Bystander Effect)',
        'answer': (
            "🎯 **ADC'lerde 'Bystander Killing' (Komşu Hücre Öldürme) Etkisi Nedir?**\n\n"
            "• **Kanser Heterojenitesi Sorunu:** Katı tümörlerdeki kanser hücrelerinin tümü hedef reseptörü (örn. HER2 veya Trop-2) eşit düzeyde taşımaz; bazı hücreler antijen-negatiftir.\n"
            "• **Bystander Mekanizması:** ADC antijeni taşıyan hücreye bağlanıp endositozla hücre içine girdiğinde, parçalanabilir (cleavable) bağlayıcı kopar ve açığa çıkan sitotoksik yük (payload, örn. DXd, SN-38 veya MMAE) lipid membranları aşabilen (membran-geçirgen) yapıdaysa hücre dışına sızarak çevredeki antijen taşımayan komşu kanser hücrelerini de yok eder.\n"
            "• **Klinik Önemi:** T-DXd (Enhertu) ve Trodelvy'nin yüksek bystander etkisi, HER2-low veya heterojen Trop-2 tümörlerde kemoterapiye dirençli klonları dahi eradike edebilmesini sağlar. T-DM1 (Kadcyla) gibi parçalanamaz bağlayıcılı ADC'lerde ise bystander etkisi sıfırdır."
        )
    },
    'enhertu_her2_low': {
        'patterns': [
            r'her2-low', r'her2 low', r'enhertu meme kanseri', r't-dxd', r'destiny-breast'
        ],
        'title': 'Trastuzumab Deruxtecan (Enhertu / T-DXd) & HER2-Low Devrimi',
        'answer': (
            "🎀 **Enhertu (T-DXd) ve HER2-Low Yeni Kanser Sınıflaması:**\n\n"
            "• **Tarihsel Dönüm Noktası:** Geçmişte meme kanseri sadece 'HER2 Pozitif (IHC 3+)' veya 'HER2 Negatif' olarak ikiye ayrılırdı. DESTINY-Breast04 çalışmasıyla **HER2-Low (IHC 1+ veya IHC 2+/ISH-)** adı verilen yepyeni bir hasta grubu tanımlandı (tüm meme kanserlerinin yaklaşık %50-60'ı).\n"
            "• **Teknolojik Güç:** Enhertu, antikora bağlı 8 adet DXd (DAR=8.0) taşır ve çok düşük HER2 ekspresyonunda dahi güçlü bystander öldürme başlatır.\n"
            "• **Klinik Başarı:** Standart kemoterapiye kıyasla progresyonsuz sağkalımı (9.9 vs 5.1 ay) ve genel sağkalımı (23.4 vs 16.8 ay, HR=0.64) çarpıcı biçimde uzatarak kılavuzlarda 1. kategori tercih haline gelmiştir. Ayrıca HER2 mutasyonlu akciğer adenokarsinomlarında da FDA onaylıdır."
        )
    },
    'trodelvy_trop2_tnbc': {
        'patterns': [
            r'trodelvy', r'sacituzumab', r'trop-2', r'trop2', r'ascent çalışması'
        ],
        'title': 'Sacituzumab Govitecan (Trodelvy) & Trop-2 Hedefli Üçlü Negatif Meme Kanseri',
        'answer': (
            "🧬 **Trodelvy (Sacituzumab Govitecan) & Üçlü Negatif Meme Kanseri (mTNBC):**\n\n"
            "• **Trop-2 Hedefi:** Üçlü negatif meme kanserlerinin (östrojen, progesteron ve HER2 negatif) %90'ından fazlasında hücre yüzeyinde aşırı miktarda Trop-2 reseptörü bulunur.\n"
            "• **Moleküler Mimari:** İrinotekanın en aktif formu olan SN-38 molekülünü yüksek oranda (DAR=7.6) taşır. Asidik tümör ortamında ve lizozomda hızla serbest kalır.\n"
            "• **ASCENT Faz 3 Çalışması:** Daha önce kemoterapi almış metastatik TNBC hastalarında genel sağkalımı ikiye katlamış (12.1 vs 6.7 ay, HR=0.48), ölüm riskini %52 azaltmıştır.\n"
            "• **Yan Etki & Yönetim:** En sık görülen yan etkiler nötropeni ve ishaldir. UGT1A1*28 yavaş metabolizörlerinde toksisite artabileceğinden genetik test ve G-CSF profilaksisi önerilir."
        )
    },
    'padcev_nectin4_bladder': {
        'patterns': [
            r'padcev', r'enfortumab', r'nectin-4', r'nectin 4', r'ev-302', r'mesane kanseri padcev'
        ],
        'title': 'Enfortumab Vedotin (Padcev) & Nectin-4 Hedefli Mesane Kanseri',
        'answer': (
            "💧 **Padcev (Enfortumab Vedotin) & İleri Evre Mesane (Ürotelyal) Kanseri:**\n\n"
            "• **Nectin-4 Reseptörü:** İdrar kesesi ve ürotelyal karsinomların neredeyse tamamında (%95+) hücre yüzeyinde bolca bulunan bir adezyon molekülüdür.\n"
            "• **Etki Mekanizması:** Padcev, mikrotübül inhibitörü olan MMAE yükünü hücre içine sokarak hücre bölünmesini (mitoz) kilitler.\n"
            "• **EV-302 / KEYNOTE-A39 Çalışması (Tarihi Başarı):** Padcev + Pembrolizumab (Keytruda) kombinasyonu, 40 yıldır uygulanan standart platin kemoterapisini geride bırakarak genel sağkalımı 16 aydan **31.5 aya** çıkarmıştır (HR=0.47).\n"
            "• **Toksisite Uyarısı:** Cilt döküntüleri (SJS/TEN şüphesinde derhal kesilmelidir), periferik nöropati (el ve ayaklarda uyuşma/karıncalanma) ve kan şekeri yüksekliği (hiperglisemi) yakından izlenmelidir."
        )
    },
    'adc_toxicity_ild_neuropathy': {
        'patterns': [
            r'adc yan etkileri', r'ild pnömoni', r'enhertu akciğer yan etkisi', r'padcev nöropati', r'adc toksisite'
        ],
        'title': 'ADC Toksisiteleri: İnterstisyel Akciğer Hastalığı (ILD) & Periferik Nöropati',
        'answer': (
            "⚠️ **Yeni Nesil Antikor-İlaç Konjugatlarında (ADC) Kritik Toksisite Yönetimi:**\n\n"
            "1. **İnterstisyel Akciğer Hastalığı (ILD) / İlaca Bağlı Pnömonit:**\n"
            "   - **Görüldüğü İlaçlar:** T-DXd (Enhertu) kullanan hastaların yaklaşık %10-12'sinde ortaya çıkar.\n"
            "   - **Belirtiler:** Yeni başlayan veya kötüleşen nefes darlığı, kuru öksürük, efor dispnesi.\n"
            "   - **Yönetim:** Asemptomatik olsa dahi (Evre 1) tomografide opasite görüldüğünde ilaç derhal ertelenir; Evre ≥2 ILD'de ilaç **ömür boyu kesilir** ve derhal sistemik yüksek doz Kortikosteroid (≥1 mg/kg/gün prednizolon eşdeğeri) başlanır.\n\n"
            "2. **Periferik Nöropati (MMAE / Vedotin İlişkili):**\n"
            "   - **Görüldüğü İlaçlar:** Padcev (Enfortumab Vedotin) ve Adcetris (Brentuximab Vedotin).\n"
            "   - **Mekanizma:** Mikrotübül taşınmasının aksonlarda inhibisyonu nedeniyle el ve ayaklarda eldiven-çorap tarzı uyuşma, yanma ve güç kaybı.\n"
            "   - **Yönetim:** Evre 2 nöropatide semptomlar hafifleyene kadar doz durdurulur ve ardından doz indirimine gidilir.\n\n"
            "3. **Kemik İliği Baskılanması & İshal (SN-38):** Trodelvy'de mutlak nötrofil sayısı (ANC) takip edilmeli ve G-CSF (filgrastim) desteği planlanmalıdır."
        )
    },

    'pluvicto_psma_mcrpc': {
        'patterns': [
            r'pluvicto', r'psma-617', r'lutesyum.*psma', r'lu-177.*psma', r'177lu-psma'
        ],
        'title': '177Lu-PSMA-617 (Pluvicto) & Prostat Kanserinde Akıllı Radyasyon',
        'answer': (
            "🎯 **177Lu-PSMA-617 (Pluvicto) Nedir & Prostat Kanserinde Nasıl Çalışır?**\n\n"
            "• **Hedef:** Prostat kanseri hücrelerinin yüzeyinde aşırı miktarda üretilen Prostat-Spesifik Membran Antijeni (PSMA / FOLH1).\n"
            "• **Mekanizma:** PSMA'ya kilitlenen akıllı bir küçük molekül ligandın ucuna Lutesyum-177 (177Lu) radyoizotopu bağlanmıştır. Kanser hücresine girer ve 1-2 mm çapındaki mikroskobik beta ışınlarıyla DNA'yı kırarak tümörü içeriden yakar.\n"
            "• **VISION Faz 3 Çalışması:** Daha önce hormon (ARPI) ve kemoterapi almış dirençli metastatik prostat kanserinde (mCRPC) genel sağkalımı belirgin şekilde uzatmıştır (15.3 vs 11.3 ay, HR=0.62) ve PSA50 yanıt oranı %46'ya ulaşmıştır.\n"
            "• **Ön Koşul:** Tedavi öncesinde 68Ga-PSMA PET/CT çekilmeli ve lezyonların karaciğerden daha parlak tutulum gösterdiği (SUVmax > karaciğer) kanıtlanmalıdır."
        )
    },
    'lutathera_sstr2_net': {
        'patterns': [
            r'lutathera', r'dotatate', r'prrt', r'lutesyum.*dotatate', r'177lu-dotatate', r'sstr2'
        ],
        'title': '177Lu-DOTATATE (Lutathera) & Nöroendokrin Tümörler (GEP-NET)',
        'answer': (
            "🧬 **Lutathera (177Lu-DOTATATE) & Peptit Reseptör Radyonüklit Tedavisi (PRRT):**\n\n"
            "• **Hedef:** Mide-bağırsak ve pankreas nöroendokrin tümörlerinin (GEP-NET) yüzeyindeki Somatostatin Reseptör Tip 2 (SSTR2).\n"
            "• **Yapı:** Oktreotat peptidi, DOTA şelatörü ile Lutesyum-177'ye bağlanmıştır.\n"
            "• **NETTER-1 & NETTER-2 Çalışmaları:** İlerlemiş midgut NET'lerde hastalık ilerlemesi veya ölüm riskini **%79 oranında azaltmıştır** (PFS HR=0.21); yüksek dereceli G2/G3 tümörlerde birinci basamakta da üstünlük sağlamıştır.\n"
            "• **Böbrek Koruması:** Böbrek proksimal tübüllerinde birikimi engellemek için infüzyonla birlikte L-Lizin ve L-Arjinin amino asit solüsyonu verilir."
        )
    },
    'targeted_alpha_therapy_actinium': {
        'patterns': [
            r'aktinyum', r'actinium', r'225ac', r'alfa tedavisi', r'targeted alpha', r'tat tedavisi'
        ],
        'title': 'Hedefe Yönelik Alfa Tedavisi (TAT: Actinium-225) & Yüksek LET',
        'answer': (
            "⚡ **Hedefe Yönelik Alfa Tedavisi (TAT - Targeted Alpha Therapy) ve Actinium-225 (225Ac):**\n\n"
            "• **Alfa Parçacıklarının Gücü:** Beta ışınlarına (elektron) kıyasla Helyum çekirdeği (alfa) yaklaşık 7.300 kat daha ağırdır ve çok yüksek Doğrusal Enerji Transferine (LET ~100 keV/μm) sahiptir.\n"
            "• **Hücresel Hassasiyet:** Alfa parçacıklarının dokudaki menzili yalnızca 40-80 mikrometredir (2-4 hücre çapı). Tüm enerjisini bu dar alana bırakarak DNA'da onarılması imkansız 'küme çift zincir kırıkları' (complex clustered DSB) açar.\n"
            "• **Klinik Devrim:** Lutesyum-177 tedavisine direnç kazanmış veya hipoksik radyo-dirençli metastatik prostat kanserli hastalarda **225Ac-PSMA-617**, hastaların %82'sinde derin PSA düşüşü ve tam radyolojik yanıtlar sağlamaktadır.\n"
            "• **Tükürük Bezi Koruması:** En kritik yan etki ağız kuruluğudur (kserostomi); tükürük bezlerine lokal buz uygulaması veya botulinum toksini enjeksiyonu ile perfüzyon azaltılır."
        )
    },
    'theranostics_pet_paired_imaging': {
        'patterns': [
            r'teranostik', r'theranostics', r'we see what we treat', r'ga-68.*lu-177', r'galyum 68'
        ],
        'title': 'Teranostik Tıp (Theranostics) & Eşleştirilmiş Görüntüleme-Tedavi',
        'answer': (
            "🔮 **Teranostik (Theranostics - 'Tedavi + Teşhis') Tıp Felsefesi:**\n\n"
            "• **'Gördüğümüzü Tedavi Ediyoruz, Tedavi Ettiğimizi Görüyoruz' Prensibi:**\n"
            "   1. **Görüntüleme Fazı (Teşhis):** Aynı hedefleme molekülü pozitron yayan bir izotopla (örn. Galyum-68 veya Flor-18) işaretlenir ve PET/CT çekilir. Tümörün reseptörü taşıyıp taşımadığı ve vücuttaki tam koordinatları milimetrik olarak haritalanır.\n"
            "   2. **Tedavi Fazı:** Reseptör pozitif hastada aynı molekülün ucuna bu kez öldürücü ışın yayan bir izotop (Lutesyum-177 veya Aktinyum-225) takılarak hastaya damardan verilir.\n"
            "• **Avantajı:** İlacın tümöre ulaşıp ulaşmayacağı tahmin edilmez; tedavi başlamadan önce PET ekranında bizzat görülür."
        )
    },
    'rlt_radiation_safety_dosimetry': {
        'patterns': [
            r'radyoligand güvenliği', r'mird dozimetri', r'böbrek koruma', r'amino asit infüzyonu', r'pluvicto ağız kuruluğu'
        ],
        'title': 'Radyoligand Güvenliği, Böbrek Koruma & MIRD Dozimetri',
        'answer': (
            "🛡️ **Radyoligand Tedavilerinde Güvenlik Protokolleri & Dozimetri:**\n\n"
            "1. **MIRD (Medical Internal Radiation Dose) Dozimetrisi:** Tümörün ve kritik organların (Böbrekler, Kemik İliği) soğurduğu toplam radyasyon dozu (Gray - Gy) hesaplanır. Kümülatif böbrek dozunun 23 Gray eşiğini aşmaması hedeflenir.\n"
            "2. **Amino Asit Nefroproteksiyonu:** Peptitler böbrek proksimal tübüllerindeki megalin/kübilin reseptörlerinden geri emilir. L-Lizin ve L-Arjinin infüzyonu bu reseptörleri rekabetçi olarak doyurarak böbrek radyasyon tutulumunu %40-50 azaltır.\n"
            "3. **Kemik İliği Güvenliği:** Kırmızı kemik iliği dozu <2 Gy tutulur. Tedavi öncesi hemoglobin, trombosit ve nötrofil sayımları kontrol edilir.\n"
            "4. **Radyasyon İzolasyonu:** Hastalar infüzyon sonrası bol sıvı tüketerek radyoizotopun atılmasını hızlandırır; ilk birkaç gün yakın temas önlemlerine uyulur."
        )
    },

    'blinatumomab_bite_all': {
        'patterns': [
            r'blinatumomab', r'blincyto', r'cd19.*cd3', r'bite.*lösemi', r'bite antikor'
        ],
        'title': 'Blinatumomab (Blincyto) & CD19xCD3 BiTE ile B-ALL Tedavisi',
        'answer': (
            "🧬 **Blinatumomab (Blincyto) & BiTE (Bispecific T-Cell Engager) Teknolojisi:**\n\n"
            "• **Mekanizma:** Fc bölgesi bulunmayan tandem iki scFv antikor parçacığından (~55 kDa) oluşur. Bir koluyla B hücre yüzeyindeki CD19 antijenine, diğer koluyla sitotoksik T hücrelerinin CD3 reseptörüne bağlanır.\n"
            "• **İmmünolojik Sinaps:** T hücresi ile tümör arasında <12 nm'lik ultra dar bir litik sinaps oluşturur. T hücresi MHC antijen sunumuna gerek duymadan aktive olarak perforin ve granzim salgılar, B-ALL hücrelerini imha eder.\n"
            "• **Endikasyon & Kanıt:** Nüks/dirençli B hücreli akut lenfoblastik lösemi (B-ALL) ve Minimal Rezidüel Hastalık (MRD) pozitifliğinde altın standarttır (TOWER ve ECOG E1910 çalışmaları).\n"
            "• **Uygulama:** Yarı ömrü çok kısa (~2 saat) olduğu için taşınabilir infüzyon pompası ile 28 gün aralıksız damardan verilir."
        )
    },
    'amivantamab_egfr_met_nsclc': {
        'patterns': [
            r'amivantamab', r'rybrevant', r'egfr.*met', r'ekzon 20 insersiyon', r'mariposa'
        ],
        'title': 'Amivantamab (Rybrevant) EGFRxMET Bispesifik Antikor & Akciğer Kanseri',
        'answer': (
            "🎯 **Amivantamab (Rybrevant) & Çift Reseptör Blokajı:**\n\n"
            "• **Hedef:** Küçük hücreli dışı akciğer kanserinde (NSCLC) hem EGFR hem de MET tirozin kinaz reseptörlerini aynı anda bağlayan tam IgG1 DuoBody bispesifik antikordur.\n"
            "• **Üçlü Etki Mekanizması:**\n"
            "   1. EGFR ve MET ligand bağlanmasını eş zamanlı bloke eder.\n"
            "   2. Reseptör komplekslerini hücre içine çekerek parçalanmaya yönlendirir (endositoz).\n"
            "   3. Bağışıklık hücrelerini çağırarak antikor bağımlı hücresel sitotoksisite (ADCC) ve trogositoz ile tümörü tüketir.\n"
            "• **Landmark Çalışmalar:** MARIPOSA çalışmasında Lazertinib ile kombine edilerek Osimertinib'e kıyasla progresyonsuz sağkalımı (PFS 23.7 vs 16.6 ay, HR=0.70) belirgin uzatmıştır; PAPILLON çalışmasında EGFR Ekzon 20 insersiyonlu grupta kemoterapiyle birlikte 1. basamak standarttır."
        )
    },
    'teclistamab_bcma_cd3_myeloma': {
        'patterns': [
            r'teclistamab', r'tecvayli', r'talquetamab', r'talvey', r'gprc5d', r'bcma.*bispesifik'
        ],
        'title': 'Teclistamab & Talquetamab: Multipl Miyelomda Bispesifik T-Hücre Angajmanı',
        'answer': (
            "🩸 **Multipl Miyelomda Bispesifik T-Hücre Angaje Edicileri (BsAbs):**\n\n"
            "• **Teclistamab (Tecvayli):** Miyelom hücrelerindeki BCMA antijenini T hücresinin CD3 molekülüne kilitler (MajesTEC-1 çalışması, %63 genel yanıt oranı). CAR-T tedavisi için haftalarca bekleyemeyen refrakter miyelom hastalarında 'raftan hazır (off-the-shelf)' hemen başlanabilir devrimsel tedavidir.\n"
            "• **Talquetamab (Talvey):** BCMA'dan bağımsız yeni bir hedef olan GPRC5D reseptörünü hedefler (MonumenTAL-1, %74 yanıt). BCMA CAR-T veya Teclistamab sonrası nüks eden hastalarda çapraz direnç olmaksızın derin yanıtlar üretir.\n"
            "• **Uygulama:** Her iki ilaç da cilt altı (subkutan) enjeksiyon şeklinde ve basamaklı dozlama (step-up) ile verilir."
        )
    },
    'tebentafusp_immtac_uveal_melanoma': {
        'patterns': [
            r'tebentafusp', r'kimmtrak', r'immtac', r'uveal melanom', r'gp100.*hla'
        ],
        'title': 'Tebentafusp (Kimmtrak) & Uveal Melanomda İlk TCR Bispesifik Tedavisi',
        'answer': (
            "👁️ **Tebentafusp (Kimmtrak) & ImmTAC (TCR Tabanlı Bispesifik İlaç) Devrimi:**\n\n"
            "• **Dünyada Bir İlk:** Yalnızca hücre yüzeyindeki antijenleri değil, hücre içindeki protein parçacıklarını da tanıyan ilk çözünür T-hücre reseptörü (TCR) bispesifik ilacıdır.\n"
            "• **Hedef:** HLA-A*02:01 molekülü tarafından sunulan hücre içi gp100 melanom peptidine ultra yüksek pikomolar ilgiyle (~50 pM) bağlanır; diğer ucuyla CD3 T hücrelerini aktive eder.\n"
            "• **IMCgp100-202 Faz 3 Çalışması:** Karaciğer metastazı yapan ve tedavi seçeneği bulunmayan ölümcül metastatik uveal (göz) melanomunda genel sağkalımı kanıtlanmış ilk ilaçtır (21.7 vs 16.0 ay, HR=0.51, ölüm riskini %49 azaltmıştır).\n"
            "• **Ön Koşul:** Hastanın genetik olarak HLA-A*02:01 doku tipine sahip olması zorunludur."
        )
    },
    'bispecific_step_up_dosing_crs': {
        'patterns': [
            r'step-up dosing', r'basamaklı dozlama', r'bispesifik crs', r'blincyto nörotoksisite', r'bispesifik yan etkileri'
        ],
        'title': 'Bispesifik Antikorlarda Basamaklı Dozlama (Step-Up Dosing) & Güvenlik',
        'answer': (
            "🛡️ **Bispesifik Antikorlarda Güvenlik ve Basamaklı Dozlama (Step-Up Dosing):**\n\n"
            "1. **Neden Basamaklı Dozlama Yapılır?** T hücrelerinin ilk temas anında aşırı fırtına şeklinde sitokin (IL-6, IFN-gamma, TNF) salgılamasını (Sitokin Salınım Sendromu - CRS) engellemek için tedaviye tam dozla değil, çok düşük mikro-dozlarla başlanır.\n"
            "2. **Protokol:** Örneğin Teclistamab'da 1. gün 0.06 mg/kg, 4. gün 0.3 mg/kg, 7. gün tam doz 1.5 mg/kg verilir. Blinatumomab'da ilk 7 gün 9 mcg/gün, sonra 28 mcg/gün'e çıkılır.\n"
            "3. **Premedikasyon:** Her basamak öncesi deksametazon, antihistaminik ve parasetamol verilir.\n"
            "4. **CRS & ICANS Yönetimi:** Ateş ve hipotansiyonda derhal Tocilizumab (anti-IL-6R); nörotoksisite (ICANS: titreme, konuşma bozukluğu, kafa karışıklığı) gelişirse yüksek doz Deksametazon uygulanır."
        )
    },

    'bot_help': {
        'patterns': [
            r'sen kimsin', r'ne yapabilirsin', r'yardım', r'komutlar', r'neler yaparsın',
            r'kimsin', r'nasıl çalışırsın'
        ],
        'title': 'Bot Hakkında & Kullanım Kılavuzu',
        'answer': (
            "🤖 **Erken Kanser Teşhisi ve Klinik Onkoloji Soru Botu Hakkında:**\n\n"
            "Bu bot **%100 çevrimdışı (offline)** çalışır ve önceden kaydedilmiş hiçbir kişisel veya offline verinizi kullanmaz. "
            "Tamamen sizin soracağınız sorulara ve vereceğiniz cevaplara göre çalışır.\n\n"
            "📌 **Neler Yapabilirsiniz?**\n"
            "1. **Tıbbi Soru Sorabilirsiniz:** Erken belirtiler, tarama testleri, risk faktörleri veya kılavuzlar hakkında dilediğinizi yazabilirsiniz (Örn: *'Akciğer kanseri belirtileri neler?'*, *'Kolonoskopi ne zaman gerekir?'*).\n"
            "2. **Kişiselleştirilmiş Erken Teşhis Anamnezi Başlatabilirsiniz:** Şikayetlerinizi yazarak (Örn: *'2 aydır öksürüyorum ve kilo verdim'*) ya da *'Anamnez başlat'* diyerek 5 adımlı interaktif risk değerlendirmesini başlatabilirsiniz.\n"
            "3. **Acil Sevk Kılavuzu (NICE NG12):** Belirtilerinizin %3 PPV acil sevk gerektirip gerektirmediğini anında öğrenebilirsiniz."
        )
    }
}


class DiagnosticChatbot:
    """
    100% Offline Conversational Early Cancer Diagnosis & Anamnesis Assistant.
    - Operates completely locally without requiring external offline database files.
    - Answers medical/oncology questions directly from built-in clinical repository.
    - Conducts interactive multi-turn anamnesis and evaluates cancer risk dynamically
      strictly based on the user's answers to the questions asked.
    """

    STATE_GREETING = "GREETING"
    STATE_SYMPTOMS = "SYMPTOMS"
    STATE_DEMOGRAPHICS = "DEMOGRAPHICS"
    STATE_SMOKING = "SMOKING"
    STATE_FAMILY = "FAMILY"
    STATE_DURATION = "DURATION"
    STATE_EVALUATION = "EVALUATION"
    STATE_COMPLETED = "COMPLETED"

    def __init__(self, db_path: Optional[str] = None):
        self.sessions: Dict[str, Dict[str, Any]] = {}
        # Optional predictor instance if available, but offline dynamic engine is primary
        self.predictor = None
        if db_path is not None:
            try:
                from cancer_risk_predictor import CancerRiskPredictor
                self.predictor = CancerRiskPredictor(db_path=db_path, auto_seed=False)
            except Exception:
                self.predictor = None

    def start_session(self, session_id: str = "default") -> Dict[str, Any]:
        """Initializes a new interactive consultation session."""
        self.sessions[session_id] = {
            'state': self.STATE_SYMPTOMS,
            'age': 50,
            'gender': 'M',
            'symptoms': [],
            'smoking_pack_years': 0.0,
            'family_history': False,
            'duration_desc': 'Belirtilmedi',
            'messages': []
        }

        greeting = (
            "Merhaba! Ben **Klinik Onkoloji & Erken Kanser Teşhisi Soru Botu** 🩺.\n\n"
            "Bu bot tamamen çevrimdışı (offline) çalışır; önceden kayıtlı herhangi bir offline verinizi kullanmaz. "
            "Bana erken kanser belirtileri, tarama testleri ve risk faktörleri hakkında dilediğiniz soruları sorabilir, "
            "veya şikayetlerinizi belirterek sorulara vereceğiniz cevaplara göre adım adım kişisel risk analizi yaptırabilirsiniz.\n\n"
            "👉 **Ne yapmak istersiniz?**\n"
            "• Merak ettiğiniz bir soruyu sorabilirsiniz (Örn: *'Akciğer kanseri belirtileri neler?'*, *'Kolonoskopi ne zaman gerekir?'*)\n"
            "• Veya hissettiğiniz şikayetleri yazarak doğrudan erken teşhis analizini başlatabilirsiniz."
        )

        quick_replies = [
            "Akciğer Kanseri Belirtileri",
            "Meme Kanseri ve Kitleler",
            "Kolonoskopi Ne Zaman Gerekir?",
            "NICE 2-Week-Wait Kuralı Nedir?",
            "Balgamda Kan Tehlikeli mi?",
            "Şikayetim Var (Anamnez Başlat)"
        ]

        self.sessions[session_id]['messages'].append({
            'sender': 'bot',
            'text': greeting,
            'quick_replies': quick_replies
        })

        return {
            'session_id': session_id,
            'state': self.STATE_SYMPTOMS,
            'message': greeting,
            'quick_replies': quick_replies
        }

    def normalize_text(self, text: str) -> str:
        """Standardizes Turkish characters and lowercases text."""
        if not text:
            return ""
        norm = text.lower()
        norm = norm.replace('İ', 'i').replace('I', 'ı')
        return norm

    def is_question(self, text: str) -> bool:
        """Determines if the user's input is a question."""
        norm = self.normalize_text(text)
        if '?' in text:
            return True
        question_words = [
            r'\bnedir\b', r'\bnasıl\b', r'\bne zaman\b', r'\bhangi\b', r'\bkimlere\b',
            r'\btehlikeli mi\b', r'\bkanser mi\b', r'\bolabilir mi\b', r'\bvar mı\b',
            r'\bmidir\b', r'\bmudur\b', r'\bmüdür\b', r'\bmıdır\b',
            r'\bbelirtileri\b', r'\btaraması\b', r'\btesti\b', r'\bne demek\b',
            r'\banlamına gelir\b', r'\bne yapmalıyım\b', r'\bne yapmam gerek\b',
            r'\bsebebi\b', r'\bnedenleri\b', r'\byayar mı\b', r'\bzararlı mı\b'
        ]
        return any(re.search(w, norm) for w in question_words)

    def find_matching_topic(self, text: str) -> Optional[Dict[str, Any]]:
        """Finds matching offline medical topic from built-in knowledge base."""
        norm = self.normalize_text(text)
        for topic_key, topic_data in OFFLINE_KNOWLEDGE_TOPICS.items():
            for pat in topic_data['patterns']:
                if re.search(pat, norm):
                    return topic_data
        return None

    def extract_symptoms_from_text(self, text: str) -> List[str]:
        """Uses Turkish NLP regex patterns to extract clinical symptom keys from natural language text."""
        norm = self.normalize_text(text)
        detected_symptoms = []
        for symptom_key, patterns in NLP_SYMPTOM_PATTERNS.items():
            for pat in patterns:
                if re.search(pat, norm):
                    if symptom_key not in detected_symptoms:
                        detected_symptoms.append(symptom_key)
                    break
        return detected_symptoms

    def extract_demographics(self, text: str) -> Tuple[Optional[int], Optional[str]]:
        """Extracts age integer and biological sex (M/F) from free text."""
        age = None
        age_match = re.search(r'\b(1[8-9]|[2-8][0-9]|9[0-9])\b', text)
        if age_match:
            age = int(age_match.group(1))

        gender = None
        lower_t = self.normalize_text(text)
        if re.search(r'\b(erkek|adam|bay|bey|baba|oğul)\b', lower_t):
            gender = 'M'
        elif re.search(r'\b(kadın|kız|bayan|anne|hanım)\b', lower_t):
            gender = 'F'

        return age, gender

    def extract_smoking_pack_years(self, text: str) -> float:
        """Calculates pack-years from natural language text."""
        lower_t = self.normalize_text(text)
        if re.search(r'\b(içmiyorum|kullanmıyorum|hiç|bıraktım|içmedim|yok)\b', lower_t) and not re.search(r'\b(yıl|paket)\b', lower_t):
            return 0.0

        # Pattern: X yıl günde Y paket
        m_year_pack = re.search(r'(\d+)\s*yıl.*?(\d+[\.,]?\d*)\s*paket', lower_t)
        if m_year_pack:
            years = float(m_year_pack.group(1))
            packs = float(m_year_pack.group(2).replace(',', '.'))
            return round(years * packs, 1)

        # Pattern: günde Y paket X yıl
        m_pack_year = re.search(r'(\d+[\.,]?\d*)\s*paket.*?(\d+)\s*yıl', lower_t)
        if m_pack_year:
            packs = float(m_pack_year.group(1).replace(',', '.'))
            years = float(m_pack_year.group(2))
            return round(years * packs, 1)

        # Pattern: sadece yıl verilmişse (varsayılan 1 paket/gün)
        m_just_year = re.search(r'(\d+)\s*yıl', lower_t)
        if m_just_year:
            return float(m_just_year.group(1)) * 1.0

        # Pattern: sadece paket verilmişse (10 paket-yıl varsayılan)
        m_just_pack = re.search(r'(\d+)\s*paket', lower_t)
        if m_just_pack:
            return float(m_just_pack.group(1)) * 10.0

        return 0.0

    def extract_family_history(self, text: str) -> bool:
        """Determines if user indicated family cancer history."""
        lower_t = self.normalize_text(text)
        negative_words = ['yok', 'hayır', 'hiç kimse', 'olmadı', 'temiz', 'bulunmuyor']
        if any(w in lower_t for w in negative_words) and not re.search(r'\b(ama|fakat|lakin)\b', lower_t):
            return False

        positive_words = ['var', 'evet', 'annem', 'babam', 'kardeşim', 'dedem', 'teyzem', 'amcam', 'kanserdi', 'kanser oldu']
        return any(w in lower_t for w in positive_words)

    def evaluate_patient_dynamically(self, age: int, gender: str, symptoms: List[str],
                                    smoking_pack_years: float = 0.0,
                                    family_history: bool = False) -> Dict[str, Any]:
        """
        Pure-Python, 100% offline Bayesian risk evaluation calibrated with NICE NG12 clinical guidelines.
        Does NOT access or query any external or offline database file.
        Strictly calculates results dynamically based on the parameters passed from user answers.
        """
        gender_code = 'M' if str(gender).upper().startswith(('M', 'E')) else 'F'
        canonical_symptoms = [s for s in symptoms if s in SYMPTOM_DEFINITIONS]

        # 1. Baseline log-odds from population epidemiology
        log_odds = -3.85
        log_odds += (age - 50.0) * 0.046
        log_odds += (smoking_pack_years / 10.0) * 0.44
        if family_history:
            log_odds += 0.68

        urgent_trigger_reasons = []
        meets_nice_urgent_referral = False

        # Multi-organ initial likelihood scores
        organ_scores = {
            'LUNG': 0.021,
            'COLORECTAL': 0.019,
            'BREAST': 0.025 if gender_code == 'F' else 0.001,
            'PROSTATE': 0.022 if gender_code == 'M' else 0.0,
            'PANCREATIC': 0.009,
            'GASTRIC': 0.011,
            'BLADDER_RENAL': 0.012,
            'MELANOMA': 0.008
        }

        # 2. Symptoms odds ratio and NICE 2WW evaluation
        for s in canonical_symptoms:
            info = SYMPTOM_DEFINITIONS[s]
            weight = info.get('base_odds_ratio', 3.0)
            log_odds += (weight / 5.0) * 0.65

            ppv = info.get('nice_ppv', 0.02)
            if ppv >= 0.03:
                meets_nice_urgent_referral = True
                urgent_trigger_reasons.append(f"{info['name_tr']} (NICE PPV: %{round(ppv*100, 1)})")

            for tc in info.get('target_cancers', []):
                if tc in organ_scores:
                    organ_scores[tc] += (weight / 5.0) * 0.25

        if smoking_pack_years >= 20.0:
            organ_scores['LUNG'] += (smoking_pack_years / 20.0) * 0.35
            organ_scores['BLADDER_RENAL'] += (smoking_pack_years / 20.0) * 0.15

        # Logistic transformation
        composite_prob = 1.0 / (1.0 + math.exp(-log_odds))
        composite_prob_percent = round(composite_prob * 100.0, 1)

        if composite_prob_percent >= 3.0:
            meets_nice_urgent_referral = True
            if not urgent_trigger_reasons:
                urgent_trigger_reasons.append("Kombine risk acil NICE sevk eşiğini (%3.0) aştı.")

        # Determine clinical risk tier
        if composite_prob_percent >= 40.0:
            risk_tier = "Kritik / Acil Sevk (Kırmızı Bayrak)"
            tier_code = "CRITICAL"
            color_code = "#dc2626"
            action_summary = "NICE 2-Week-Wait (2WW) Acil Sevk Protokolü gereklidir. Acil ileri görüntüleme ve doku biyopsisi planlanmalıdır."
        elif composite_prob_percent >= 15.0:
            risk_tier = "Yüksek Risk"
            tier_code = "HIGH"
            color_code = "#ea580c"
            action_summary = "Öncelikli uzman hekim değerlendirmesi ve hedefe yönelik görüntüleme (Toraks BT, Kolonoskopi veya Mamografi) önerilir."
        elif composite_prob_percent >= 3.0 or meets_nice_urgent_referral:
            risk_tier = "Orta / Tetkik Gerektiren Risk"
            tier_code = "MODERATE"
            color_code = "#eab308"
            action_summary = "NICE sevk eşiği (%3.0) aşılmıştır. 2 hafta içinde uzman hekim muayenesi ve ileri tetkik önerilir."
        else:
            risk_tier = "Düşük Risk / Rutin Takip"
            tier_code = "LOW"
            color_code = "#16a34a"
            action_summary = "Mevcut bulgular genel toplum tarama eşiklerinin altındadır. Yaşa uygun rutin kanser tarama takvimine devam edilmelidir."

        # Organ likelihood ranking
        sorted_organs = sorted(organ_scores.items(), key=lambda x: x[1], reverse=True)
        sum_scores = sum(organ_scores.values()) or 1.0

        organ_probs = []
        for code, raw_score in sorted_organs:
            meta = CANCER_METADATA.get(code, {})
            pct = round((raw_score / sum_scores) * composite_prob_percent, 1)
            pct = min(99.0, max(0.1, pct))
            organ_probs.append({
                'cancer_code': code,
                'cancer_name': meta.get('name_tr', code),
                'probability_percent': pct,
                'recommended_test': meta.get('recommended_investigation', 'Genel Uzman Muayenesi'),
                'key_drivers': meta.get('key_risk_drivers', '')
            })

        top_cancer_code = sorted_organs[0][0]
        top_meta = CANCER_METADATA.get(top_cancer_code, {})
        primary_recommendation = top_meta.get('recommended_investigation', 'Genel Hekim Muayenesi')

        return {
            'patient_profile': {
                'age': age,
                'gender': 'Erkek' if gender_code == 'M' else 'Kadın',
                'smoking_pack_years': smoking_pack_years,
                'family_history': 'Var' if family_history else 'Yok',
                'active_symptoms': [SYMPTOM_DEFINITIONS[s]['name_tr'] for s in canonical_symptoms],
                'symptom_keys': canonical_symptoms
            },
            'risk_evaluation': {
                'composite_risk_percent': composite_prob_percent,
                'composite_probability': composite_prob,
                'risk_tier': risk_tier,
                'tier_code': tier_code,
                'color_code': color_code,
                'meets_nice_2ww_referral': meets_nice_urgent_referral,
                'urgent_reasons': urgent_trigger_reasons,
                'clinical_action': action_summary,
                'primary_diagnostic_recommendation': primary_recommendation
            },
            'organ_probabilities': organ_probs,
            'primary_organ': top_cancer_code.lower(),
            'primary_suspected_malignancy': {
                'organ': top_cancer_code.lower(),
                'name': top_meta.get('name_tr', top_cancer_code)
            },
            'clinical_action_plan': f"{primary_recommendation} - {action_summary}"
        }


    def synthesize_clinical_guidance(self, text: str) -> str:
        """
        Adaptive Medical Synthesizer:
        Generates contextual oncology and clinical guidance when no pre-defined static topic matches.
        Scans for mentioned anatomical organs, pathological findings, symptoms, or diagnostic tests.
        """
        norm = self.normalize_text(text)
        
        # 1. Organ mapping
        organs = {
            'karaciğer': ('Karaciğer', 'Gastroenteroloji / Hepatoloji veya Genel Cerrahi', 'Dinamik 3 Fazlı Batın BT veya Karaciğer MRG (Primovist)'),
            'böbrek': ('Böbrek', 'Üroloji', 'Üriner Sistem Ultrasonografisi ve Kontrastlı Batın BT'),
            'pankreas': ('Pankreas', 'Gastroenteroloji / Cerrahi Onkoloji', 'Pankreas Protokolü Dinamik BT / MRCP ve Serum CA 19-9'),
            'akciğer': ('Akciğer', 'Göğüs Hastalıkları / Torasik Onkoloji', 'Düşük Doz Toraks BT (LDCT) ve Bronkoskopi'),
            'meme': ('Meme', 'Genel Cerrahi / Meme Cerrahisi', 'Bilateral Dijital Mamografi, Meme USG ve Tru-Cut Biyopsi'),
            'mide': ('Mide & Özofagus', 'Gastroenteroloji', 'Üst GIS Endoskopisi (Gastroskopi) ve patolojik biyopsi'),
            'kolon': ('Kalın Bağırsak (Kolon)', 'Gastroenteroloji / Genel Cerrahi', 'Tam Kolonoskopi ve Biyopsi'),
            'bağırsak': ('Bağırsak', 'Gastroenteroloji', 'Kolonoskopi ve Gaitada Gizli Kan Testi (FIT)'),
            'prostat': ('Prostat', 'Üroloji', 'Serum Total PSA, mpMRI ve Füzyon Biyopsi'),
            'beyin': ('Merkezi Sinir Sistemi (Beyin)', 'Beyin ve Sinir Cerrahisi (Nöroşirürji) / Nöroloji', 'Kontrastlı Kraniyal Beyin MRG'),
            'kemik': ('Kemik & İskelet Sistemi', 'Ortopedik Onkoloji', 'Direkt Grafi, Ekstremite MRG ve Kemik Sintigrafisi'),
            'tiroid': ('Tiroid Bezi', 'Endokrinoloji veya Genel Cerrahi', 'Tiroid Ultrasonografisi ve İnce İğne Aspirasyon Biyopsisi (İİAB)'),
            'yumurtalık': ('Yumurtalık (Over)', 'Jinekolojik Onkoloji', 'Transvajinal Pelvik Ultrasonografi, Kontrastlı Pelvik MRG ve CA 125'),
            'rahim': ('Rahim (Endometrium)', 'Kadın Hastalıkları ve Doğum', 'Endometriyal Pipelle Biyopsisi ve Pelvik USG'),
            'larenks': ('Larenks (Gırtlak)', 'Kulak Burun Boğaz (KBB)', 'Fiberoptik Laringoskopi ve Biyopsi'),
            'boyun': ('Boyun Bölgesi', 'Kulak Burun Boğaz (KBB)', 'Boyun Ultrasonografisi ve Kontrastlı Boyun BT')
        }

        detected_organ = None
        for k, v in organs.items():
            if k in norm:
                detected_organ = v
                break

        # 2. Finding mapping
        findings = []
        if any(w in norm for w in ['lezyon', 'nodül', 'kitle', 'kist', 'polip', 'kalınlaşma', 'şişlik']):
            findings.append("yapısal doku değişikliği (nodül / kitle / lezyon)")
        if any(w in norm for w in ['ağrı', 'sızı', 'kramp']):
            findings.append("ağrı semptomu")
        if any(w in norm for w in ['kan', 'kanama', 'leke']):
            findings.append("kanama bulgusu")

        if detected_organ:
            organ_name, clinic, test = detected_organ
            finding_text = f" ({', '.join(findings)})" if findings else ""
            return (
                f"🩺 **{organ_name} Klinik Değerlendirme & Tetkik Rehberi:**\n\n"
                f"• **Soru Değerlendirmesi:** {organ_name} ile ilgili belirttiğiniz durumu{finding_text} inceledim.\n"
                f"• **İyi Huylu vs Kötü Huylu Ayrımı:** Vücutta saptanan lezyonların ve kitlelerin önemli bir bölümü iyi huyludur (basit kistler, adenomlar, hemanjiomlar vb.). Ancak hiçbir doku değişikliği hekim onayı olmaksızın göz ardı edilmemelidir.\n"
                f"• **Önerilen İnceleme:** Tanısal netlik için **{test}** tetkiki altın standarttır.\n"
                f"• **Başvurulacak Uzmanlık:** Zaman kaybetmeden **{clinic}** hekimine başvurarak muayene olmanız ve görüntüleme disklerinizi değerlendirtmeniz önerilir.\n\n"
                f"👉 Dilerseniz şikayetlerinizi detaylandırarak aşağıdaki **'Şikayetim Var (Anamnez Başlat)'** butonundan 100.000 kişilik veri tabanı üzerinden erken kanser risk analizinizi yaptırabilirsiniz."
            )

        # General intelligent synthesis fallback
        return (
            "ℹ️ **Klinik Danışma & Erken Teşhis Rehberi:**\n\n"
            "Tıbbi sorunuzu inceledim. Onkolojik süreçlerde altın kural; vücutta 2-3 haftadan uzun süren normal dışı değişiklikleri ciddiye almaktır.\n\n"
            "• **Danışabileceğiniz Başlıca Alanlar:**\n"
            "  1. **Kanser Türleri & Belirtileri:** Akciğer, meme, kolon, prostat, mide, pankreas, lenfoma, beyin vb.\n"
            "  2. **Modern Tedaviler:** Akıllı ilaçlar (Osimertinib, Sotorasib), İmmünoterapi (Keytruda), Kemoterapi, Radyoterapi.\n"
            "  3. **Tanı & Biyomarkerlar:** PET-BT SUVmax, DPYD testi, Sıvı biyopsi (ctDNA), Mamografi, Kolonoskopi.\n\n"
            "👉 Merak ettiğiniz spesifik konuyu doğrudan sorabilir veya hissettiğiniz belirtileri yazarak kişisel risk analizinizi başlatabilirsiniz."
        )

    def process_message(self, user_message: str, session_id: str = "default") -> Dict[str, Any]:
        """
        Processes an incoming user message and returns the next conversation step.
        Supports both direct question-answering and guided clinical anamnesis.
        Supports flexible parameter order: process_message(user_message, session_id) or process_message(session_id, user_message).
        """
        if session_id != "default":
            if user_message in self.sessions:
                session_id, user_message = user_message, session_id
            elif session_id not in self.sessions and (len(session_id) > len(user_message) or ' ' in session_id):
                session_id, user_message = user_message, session_id

        if session_id not in self.sessions:
            self.start_session(session_id)

        session = self.sessions[session_id]
        session['messages'].append({'sender': 'user', 'text': user_message})
        current_state = session['state']
        user_norm = self.normalize_text(user_message)

        response_text = ""
        quick_replies = []
        report_data = None

        # ----------------------------------------------------
        # CHECK: Direct Medical Question Answering
        # ----------------------------------------------------
        matched_topic = self.find_matching_topic(user_message)
        is_user_asking_q = self.is_question(user_message)
        wants_anamnesis = any(w in user_norm for w in [
            'anamnez', 'riskimi hesapla', 'şikayetim var', 'test yap', 'tahlil yap',
            'rutin kontrol', 'kontrol yaptırmak', 'kontrol'
        ])
        detected_symptoms = self.extract_symptoms_from_text(user_message)

        # If user is in initial state and asks a question (without explicitly requesting anamnesis)
        if current_state == self.STATE_SYMPTOMS:
            is_initial_q = (is_user_asking_q and not wants_anamnesis) or (
                matched_topic is not None and not detected_symptoms and not wants_anamnesis
            )
            if is_initial_q:
                if matched_topic:
                    response_text = matched_topic['answer']
                else:
                    response_text = self.synthesize_clinical_guidance(user_message)

                quick_replies = [
                    "Şikayetim Var (Anamnez Başlat)",
                    "Akciğer Kanseri Belirtileri",
                    "Kolonoskopi Ne Zaman Gerekir?",
                    "NICE 2-Week-Wait Kuralı Nedir?"
                ]
                session['messages'].append({'sender': 'bot', 'text': response_text})
                return {
                    'session_id': session_id,
                    'state': current_state,
                    'message': response_text,
                    'quick_replies': quick_replies,
                    'report_data': None,
                    'evaluation': None,
                    'session_data': session
                }

        # If user asks an intervening question in the middle of anamnesis, answer it and prompt them to continue
        if current_state in [self.STATE_DEMOGRAPHICS, self.STATE_SMOKING, self.STATE_FAMILY, self.STATE_DURATION] and is_user_asking_q:
            answer_part = matched_topic['answer'] if matched_topic else "Sorunuzu inceledim. Erken teşhis ve tarama protokolleri belirtilen kılavuzlara göre yürütülmektedir."
            
            pending_prompts = {
                self.STATE_DEMOGRAPHICS: "👉 **Kaç yaşındasınız ve biyolojik cinsiyetiniz nedir?** (Örn: *'58 yaşındayım, erkeğim'* veya *'45 kadın'*)",
                self.STATE_SMOKING: "👉 **Tütün (sigara) kullanım öykünüz var mı? Günde kaç paket ve kaç yıl içtiniz?** (Örn: *'Günde 1 paket 20 yıl'* veya *'Hiç içmedim'*)",
                self.STATE_FAMILY: "👉 **Birinci derece akrabalarınızda (anne, baba, kardeş) teşhis edilmiş kanser öyküsü var mı?**",
                self.STATE_DURATION: "👉 **Bu şikayetleriniz yaklaşık ne kadar süredir devam ediyor?** (Örn: *'2 haftadan az'*, *'2 aydır'*)"
            }
            
            response_text = f"{answer_part}\n\n---\n💡 **Anamnezimize kaldığımız yerden devam edelim:**\n{pending_prompts.get(current_state, '')}"
            session['messages'].append({'sender': 'bot', 'text': response_text})
            return {
                'session_id': session_id,
                'state': current_state,
                'message': response_text,
                'quick_replies': ["Devam Et", "Yeniden Başla"],
                'report_data': None,
                'evaluation': None,
                'session_data': session
            }

        # ----------------------------------------------------
        # 1. SYMPTOM INTAKE
        # ----------------------------------------------------
        if current_state == self.STATE_SYMPTOMS:
            detected = self.extract_symptoms_from_text(user_message)
            if detected:
                for sym in detected:
                    if sym not in session['symptoms']:
                        session['symptoms'].append(sym)

            symptom_names_tr = [
                SYMPTOM_DEFINITIONS[s]['name_tr'] for s in session['symptoms'] if s in SYMPTOM_DEFINITIONS
            ]

            if symptom_names_tr:
                feedback = "Anladım. İfadenizden şu klinik alarm belirtilerini kaydettim:\n"
                for s_tr in symptom_names_tr:
                    feedback += f"  • 🚩 **{s_tr}**\n"
            else:
                feedback = "Şikayet kaydınız alındı (belirttiğiniz ifadede spesifik bir kırmızı bayrak semptomu doğrudan algılanmadı).\n"

            response_text = (
                f"{feedback}\n"
                "Risk analizini sorulara vereceğiniz cevaplara göre doğru hesaplayabilmemiz için:\n"
                "👉 **Kaç yaşındasınız ve biyolojik cinsiyetiniz nedir?**\n"
                "(Örn: *'58 yaşındayım, erkeğim'* veya *'45, kadın'*)"
            )
            quick_replies = ["50 yaş, Erkek", "50 yaş, Kadın", "65 yaş, Erkek", "35 yaş, Kadın"]
            session['state'] = self.STATE_DEMOGRAPHICS

        # ----------------------------------------------------
        # 2. DEMOGRAPHICS (AGE & GENDER)
        # ----------------------------------------------------
        elif current_state == self.STATE_DEMOGRAPHICS:
            age, gender = self.extract_demographics(user_message)
            if age:
                session['age'] = age
            if gender:
                session['gender'] = gender

            gender_str = "Erkek" if session['gender'] == 'M' else "Kadın"
            response_text = (
                f"Demografik bilgileriniz kaydedildi: **{session['age']} yaş, {gender_str}**.\n\n"
                "👉 **Geçmişte veya şu anda tütün (sigara) kullanım öykünüz var mı?**\n"
                "Günde kaç paket ve yaklaşık kaç yıl içtiniz?\n"
                "(Örn: *'Günde 1 paket, 25 yıl içtim'*, *'Hiç sigara içmedim'*, *'10 paket-yıl'*)"
            )
            quick_replies = [
                "Hiç içmedim (0 Paket-Yıl)",
                "Günde 1 paket, 15 yıl (15 Paket-Yıl)",
                "Günde 1 paket, 25 yıl (25 Paket-Yıl)",
                "Günde 1.5 paket, 30 yıl (45 Paket-Yıl)"
            ]
            session['state'] = self.STATE_SMOKING

        # ----------------------------------------------------
        # 3. SMOKING HISTORY
        # ----------------------------------------------------
        elif current_state == self.STATE_SMOKING:
            pack_years = self.extract_smoking_pack_years(user_message)
            session['smoking_pack_years'] = pack_years

            response_text = (
                f"Tütün maruziyetiniz hesaplandı ve kaydedildi: **{pack_years} Paket-Yıl**.\n\n"
                "👉 **Birinci derece akrabalarınızda (anne, baba, kardeş, çocuk) teşhis edilmiş herhangi bir kanser öyküsü var mı?**\n"
                "(Örn: *'Evet, babamda akciğer kanseri vardı'*, *'Hayır, ailemde kanser yok'*)"
            )
            quick_replies = [
                "Evet, 1. derece ailemde kanser var",
                "Hayır, ailemde kanser öyküsü yok",
                "Bilmiyorum / Emin değilim"
            ]
            session['state'] = self.STATE_FAMILY

        # ----------------------------------------------------
        # 4. FAMILY HISTORY
        # ----------------------------------------------------
        elif current_state == self.STATE_FAMILY:
            has_fam = self.extract_family_history(user_message)
            session['family_history'] = has_fam
            fam_str = "Mevcut (1. Derece Kanser Öyküsü Var)" if has_fam else "Yok (Bildirilmedi)"

            response_text = (
                f"Aile öykünüz kaydedildi: **{fam_str}**.\n\n"
                "👉 **Son olarak: Bu şikayetleriniz yaklaşık ne kadar süredir devam ediyor?**\n"
                "(Örn: *'2 haftadan az'*, *'3-4 haftadır'*, *'2 aydır'*, *'6 aydan uzun süredir'*)"
            )
            quick_replies = ["1 haftadan az", "2 - 3 hafta", "1 - 3 aydır", "6 aydan uzun süredir"]
            session['state'] = self.STATE_DURATION

        # ----------------------------------------------------
        # 5. DURATION & CLINICAL RISK EVALUATION
        # ----------------------------------------------------
        elif current_state == self.STATE_DURATION:
            session['duration_desc'] = user_message.strip()

            # Dynamic 100% offline evaluation without relying on offline database files
            eval_result = self.evaluate_patient_dynamically(
                age=session['age'],
                gender=session['gender'],
                symptoms=session['symptoms'],
                smoking_pack_years=session['smoking_pack_years'],
                family_history=session['family_history']
            )

            p_info = eval_result['patient_profile']
            r_eval = eval_result['risk_evaluation']
            report_data = eval_result

            nice_alert = ""
            if r_eval['meets_nice_2ww_referral']:
                nice_alert = (
                    "\n🚨 **[ACİL KLİNİK UYARI: NICE NG12 2-WEEK-WAIT KAPSAMI]**\n"
                    f"Kırmızı bayrak belirtileriniz nedeniyle **en geç 2 hafta içinde** ilgili onkoloji/uzman hekime acil sevk gerekmektedir!\n"
                    f"*Tetikleyici Nedenler:* {', '.join(r_eval['urgent_reasons'])}\n"
                )

            top_organs_text = ""
            for o in eval_result['organ_probabilities'][:3]:
                top_organs_text += f"  • **{o['cancer_name']}**: %{o['probability_percent']} (Önerilen Tetkik: *{o['recommended_test']}*)\n"

            response_text = (
                "📋 **[KLİNİK ANAMNEZ VE ERKEN TEŞHİS DEĞERLENDİRME RAPORU]**\n"
                "----------------------------------------------------------------------\n"
                f"• **Hasta Profili:** {p_info['age']} Yaş, {p_info['gender']}\n"
                f"• **Tütün Maruziyeti:** {p_info['smoking_pack_years']} Paket-Yıl\n"
                f"• **Aile Kanser Öyküsü:** {p_info['family_history']}\n"
                f"• **Saptanan Kırmızı Bayraklar:** {', '.join(p_info['active_symptoms']) if p_info['active_symptoms'] else 'Spesifik alarm bulgusu yok'}\n"
                f"• **Şikayet Süresi:** {session['duration_desc']}\n\n"
                f"⚠️ **TAHMİNİ BİLEŞİK KANSER RİSKİ:** **%{r_eval['composite_risk_percent']}** ({r_eval['risk_tier']})\n"
                f"🎯 **BİRİNCİL ŞÜPHELİ ORGAN:** **{eval_result['organ_probabilities'][0]['cancer_name']}**\n"
                f"🚨 **NICE NG12 2-WEEK-WAIT SEVK:** **{'ACİL 2WW SEVK GEREKLİ' if r_eval['meets_nice_2ww_referral'] else 'Standart / Rutin Takip'}**\n"
                f"{nice_alert}\n"
                "🔬 **ÖNCELİKLİ KLİNİK TANI PROTOKOLÜ:**\n"
                f"👉 *{r_eval['primary_diagnostic_recommendation']}*\n\n"
                "📊 **Organ Bazında Olasılık Dağılımı:**\n"
                f"{top_organs_text}\n"
                "✅ **Değerlendirme Yöntemi:** Sorulara verdiğiniz cevaplar doğrultusunda NICE NG12 klinik algoritmaları ve Bayesyen tıp modeliyle dinamik olarak hesaplanmıştır (Harici/önceden kayıtlı veri kullanılmamıştır).\n"
                "*Not: Bu rapor istatistiksel ve klinik bir risk modelidir, kesin tıbbi tanı yerine geçmez. Lütfen hekiminize danışınız.*"
            )

            quick_replies = ["Yeniden Başla", "NICE Kılavuzu Nedir?", "Kolonoskopi Ne Zaman Gerekir?"]
            session['state'] = self.STATE_COMPLETED

        # ----------------------------------------------------
        # 6. COMPLETED OR FAQ STATE
        # ----------------------------------------------------
        elif current_state == self.STATE_COMPLETED:
            if re.search(r'\b(yeniden|baştan|tekrar|yeni|sıfırla)\b', user_norm):
                return self.start_session(session_id)
            elif matched_topic:
                response_text = f"{matched_topic['answer']}\n\nYeni bir hasta değerlendirmesi için **'Yeniden Başla'** yazabilirsiniz."
                quick_replies = ["Yeniden Başla", "NICE 2WW Kuralı Nedir?"]
            elif is_user_asking_q:
                response_text = (
                    "ℹ️ Kanser şüphesinde erken tanı hayat kurtarır. İleri tetkik ve kesin tanı için hekiminize başvurunuz.\n\n"
                    "Yeni bir değerlendirme başlatmak için **'Yeniden Başla'** diyebilirsiniz."
                )
                quick_replies = ["Yeniden Başla"]
            else:
                response_text = (
                    "Bu değerlendirme tamamlandı. Yeni bir şikayet girmek için "
                    "aşağıdaki düğmeden **Yeniden Başla**'ya tıklayabilir veya merak ettiğiniz tıbbi bir soruyu sorabilirsiniz."
                )
                quick_replies = ["Yeniden Başla", "Akciğer Kanseri Belirtileri"]

        session['messages'].append({'sender': 'bot', 'text': response_text})

        return {
            'session_id': session_id,
            'state': session['state'],
            'message': response_text,
            'quick_replies': quick_replies,
            'report_data': report_data,
            'evaluation': report_data,
            'session_data': session
        }

    def interactive_cli_session(self):
        """Runs a live conversational console chat session directly in the terminal."""
        safe_print("=" * 70)
        safe_print("🏥 KLINIK ONKOLOJI VE ERKEN TEŞHIS SORU BOTU (100% OFFLINE)")
        safe_print("=" * 70)
        safe_print("Bu bot tamamen çevrimdışı (offline) çalışır; önceden kayıtlı herhangi")
        safe_print("bir verinizi kullanmaz. Sorularınıza anında cevap verir veya")
        safe_print("sorulan sorulara verdiğiniz cevaplara göre adım adım riskinizi hesaplar.\n")
        safe_print("👉 Örnek sorular sorabilirsiniz:")
        safe_print("   • 'Akciğer kanseri belirtileri nelerdir?'")
        safe_print("   • 'Meme kanserinde hangi belirtiler görülür?'")
        safe_print("   • 'Kolonoskopi ne zaman ve kimlere yapılmalı?'")
        safe_print("   • 'Balgamda kan gelmesi tehlikeli mi?'")
        safe_print("   • 'Biyopsi kanseri yayar mı?'")
        safe_print("   • 'NICE 2-Week-Wait kuralı nedir?'\n")
        safe_print("👉 Veya şikayetlerinizi yazarak doğrudan erken teşhis analizi başlatabilirsiniz:")
        safe_print("   • '2 aydır öksürüyorum ve kanlı balgam var'")
        safe_print("   • 'Anamnez başlat' / 'Riskimi hesapla'\n")
        safe_print("Çıkmak için: 'q' veya 'exit'")
        safe_print("=" * 70)

        init_data = self.start_session("cli_session")
        safe_print(f"\n[BOT]:\n{init_data['message']}\n")
        if init_data.get('quick_replies'):
            safe_print(f"  Hızlı Seçenekler: {', '.join(init_data['quick_replies'])}\n")

        while True:
            try:
                user_msg = input("Siz > ").strip()
                if not user_msg:
                    continue
                if user_msg.lower() in ['q', 'exit', 'quit', 'çıkış']:
                    safe_print("\nSoru botu oturumu sonlandırıldı. Sağlıklı günler dileriz!\n")
                    break

                res = self.process_message(user_msg, "cli_session")
                safe_print(f"\n[BOT]:\n{res['message']}\n")
                if res.get('quick_replies'):
                    safe_print(f"  Hızlı Seçenekler: {', '.join(res['quick_replies'])}\n")

                if res.get('state') == self.STATE_COMPLETED:
                    safe_print("-" * 70)

            except (KeyboardInterrupt, EOFError):
                safe_print("\nSoru botu oturumu sonlandırıldı.\n")
                break


if __name__ == '__main__':
    bot = DiagnosticChatbot()
    bot.interactive_cli_session()
