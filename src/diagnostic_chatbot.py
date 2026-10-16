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
            r'dışkı.*kan', r'büyük abdest.*kan', r'kolonoskopi', r'makat.*kan', r'rektal kanama'
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
            r'mide.*kanser', r'gastrik.*kanser', r'mide.*belirti', r'yutma güç',
            r'disfaji', r'erken doyma', r'mide.*ağrı', r'endoskopi'
        ],
        'title': 'Mide Kanseri & Disfaji (Yutma Güçlüğü)',
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
            r'mesane.*kanser', r'idrar.*kan', r'hematüri', r'böbrek.*kanser',
            r'idrarda.*kırmızı', r'mesane.*belirti', r'sistoskopi'
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
            r'melanom', r'cilt.*kanser', r'deri.*kanser', r'ben.*değiş',
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
            r'kolonoskopi.*nedir', r'kolonoskopi.*acı', r'kolonoskopi.*nasıl',
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
            r'nedensiz kilo', r'açıklanamayan kilo', r'hızlı zayıflama', r'istemsiz kilo', r'kilo verdim.*neden'
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
                    response_text = (
                        "ℹ️ **Klinik Danışma Yanıtı:**\n"
                        "Kanser erken teşhisinde altın kural; vücutta 2-3 haftadan uzun süren normal dışı değişiklikleri ciddiye almaktır.\n\n"
                        "• **En Önemli Kırmızı Bayraklar:**\n"
                        "  1. Açıklanamayan kilo kaybı (>%10) ve iştahsızlık.\n"
                        "  2. Memede, koltuk altında veya boyunda ele gelen ağrısız kitleler.\n"
                        "  3. Balgamda, idrarda veya dışkıda kan görülmesi.\n"
                        "  4. 3 haftayı geçen inatçı öksürük veya ses kısıklığı.\n"
                        "  5. Yutma güçlüğü veya bağırsak alışkanlığında kalıcı değişiklikler.\n\n"
                        "👉 Kişiselleştirilmiş bir risk değerlendirmesi yaptırmak isterseniz kendi şikayetlerinizi buraya yazabilir veya aşağıdaki **'Şikayetim Var (Anamnez Başlat)'** butonunu seçebilirsiniz."
                    )

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
