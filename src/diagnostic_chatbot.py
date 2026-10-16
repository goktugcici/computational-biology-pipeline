"""
Diagnostic Chatbot & Early Cancer Detection Anamnesis Assistant.
Provides multi-turn conversational question-and-answer symptom intake,
Turkish Natural Language Processing (NLP) symptom extraction,
NICE NG12 clinical alarm evaluation, and 100,000+ patient cohort matching.
"""

import re
import json
from typing import Dict, List, Any, Optional, Tuple
from cancer_risk_predictor import CancerRiskPredictor, SYMPTOM_DEFINITIONS

# Turkish NLP Keyword Mapping to 14 NICE NG12 Standardized Symptoms
NLP_SYMPTOM_PATTERNS: Dict[str, List[str]] = {
    'hemoptysis': [
        r'kan.*balgam', r'balgam.*kan', r'kan.*tükür', r'tükür.*kan',
        r'ağzı.*kan', r'öksür.*kan', r'hemoptizi', r'kanlı balgam'
    ],
    'chronic_cough': [
        r'öksür', r'öksürük', r'öksürüyorum', r'ses kısıklığı', r'boğazımda gıcık',
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
    ]
}


class DiagnosticChatbot:
    """
    Conversational Early Cancer Diagnosis & Anamnesis Question Bot.
    Guides the user through clinically validated question rounds,
    extracts symptoms using natural language processing,
    and returns 100k cohort evidence-based risk evaluations.
    """

    STATE_GREETING = "GREETING"
    STATE_SYMPTOMS = "SYMPTOMS"
    STATE_DEMOGRAPHICS = "DEMOGRAPHICS"
    STATE_SMOKING = "SMOKING"
    STATE_FAMILY = "FAMILY"
    STATE_DURATION = "DURATION"
    STATE_EVALUATION = "EVALUATION"
    STATE_COMPLETED = "COMPLETED"

    def __init__(self):
        self.predictor = CancerRiskPredictor()
        self.sessions: Dict[str, Dict[str, Any]] = {}

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
            "Merhaba! Ben **Erken Kanser Teşhisi ve Klinik Triyaj Soru Botu** 🩺.\n\n"
            "Amacım, hissettiğiniz belirtileri adım adım dinleyerek analiz etmek ve "
            "**100.000 kişilik doğrulanmış epidemiyoloji veri tabanımız** ile NICE NG12 klinik kılavuzları "
            "ışığında erken risk profilinizi değerlendirmektir.\n\n"
            "👉 **Şu anda sizi rahatsız eden, vücudunuzda hissettiğiniz şikayet veya belirtiler nelerdir?**\n"
            "(Örn: *'2 aydır öksürüyorum ve balgamda kan gördüm'*, *'son zamanlarda istemsiz kilo verdim'*, "
            "*'büyük abdestimde kan var'*, *'yutkunurken yemekler takılıyor'* gibi kendi cümlelerinizle yazabilirsiniz)"
        )

        self.sessions[session_id]['messages'].append({
            'sender': 'bot',
            'text': greeting,
            'quick_replies': [
                "Kanlı Balgam ve Öksürük",
                "Açıklanamayan Kilo Kaybı",
                "Dışkıda Kan / Rektal Kanama",
                "Memede Ele Gelen Kitle",
                "Yutma Güçlüğü",
                "Gözlerde/Ciltte Sarılık",
                "Hiçbir Şikayetim Yok (Rutin Kontrol)"
            ]
        })

        return {
            'session_id': session_id,
            'state': self.STATE_SYMPTOMS,
            'message': greeting,
            'quick_replies': self.sessions[session_id]['messages'][-1]['quick_replies']
        }

    def extract_symptoms_from_text(self, text: str) -> List[str]:
        """
        Uses Turkish NLP regex patterns to extract clinical symptom keys from natural language text.
        """
        normalized_text = text.lower()
        # Turkish character standardizations
        normalized_text = normalized_text.replace('İ', 'i').replace('I', 'ı')

        detected_symptoms = []
        for symptom_key, patterns in NLP_SYMPTOM_PATTERNS.items():
            for pat in patterns:
                if re.search(pat, normalized_text):
                    if symptom_key not in detected_symptoms:
                        detected_symptoms.append(symptom_key)
                    break

        return detected_symptoms

    def extract_demographics(self, text: str) -> Tuple[Optional[int], Optional[str]]:
        """Extracts age integer and gender (M/F) from free text."""
        # Age extraction
        age = None
        age_match = re.search(r'\b(1[8-9]|[2-8][0-9]|9[0-9])\b', text)
        if age_match:
            age = int(age_match.group(1))

        # Gender extraction
        gender = None
        lower_t = text.lower()
        if re.search(r'\b(erkek|adam|bay|bey|baba|oğul)\b', lower_t):
            gender = 'M'
        elif re.search(r'\b(kadın|kız|bayan|anne|hanım)\b', lower_t):
            gender = 'F'

        return age, gender

    def extract_smoking_pack_years(self, text: str) -> float:
        """Calculates pack-years from natural language text."""
        lower_t = text.lower()
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
        lower_t = text.lower()
        negative_words = ['yok', 'hayır', 'hiç kimse', 'olmadı', 'temiz', 'bulunmuyor']
        if any(w in lower_t for w in negative_words) and not re.search(r'\b(ama|fakat|lakin)\b', lower_t):
            return False

        positive_words = ['var', 'evet', 'annem', 'babam', 'kardeşim', 'dedem', 'teyzem', 'amcam', 'kanserdi', 'kanser oldu']
        return any(w in lower_t for w in positive_words)

    def process_message(self, user_message: str, session_id: str = "default") -> Dict[str, Any]:
        """
        Processes an incoming user message and returns the next conversation step.
        Supports both process_message(user_message, session_id) and process_message(session_id, user_message).
        """
        if user_message in self.sessions and session_id != "default":
            # Swapped arguments: process_message(session_id, user_message)
            session_id, user_message = user_message, session_id
        elif session_id not in self.sessions and user_message not in self.sessions:
            if session_id == "default" and "patient_" in user_message:
                pass

        if session_id not in self.sessions:
            self.start_session(session_id)

        session = self.sessions[session_id]
        session['messages'].append({'sender': 'user', 'text': user_message})
        current_state = session['state']

        response_text = ""
        quick_replies = []
        report_data = None

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
                feedback = "Anladım. İfadenizden şu klinik belirtileri kaydettim:\n"
                for s_tr in symptom_names_tr:
                    feedback += f"  • 🚩 **{s_tr}**\n"
            else:
                feedback = "Şikayet kaydınız alındı (veya şu an belirgin bir akut alarm semptomu ifade etmediniz).\n"

            response_text = (
                f"{feedback}\n"
                "Risk analizini demografik olarak doğru hesaplayabilmemiz için:\n"
                "👉 **Kaç yaşındasınız ve biyolojik cinsiyetiniz nedir?**\n"
                "(Örn: *'58 yaşındayım, erkeğim'* veya *'45, kadın'*)"
            )
            quick_replies = ["50 yaş, Erkek", "50 yaş, Kadın", "65 yaş, Erkek", "35 yaş, Kadın"]
            session['state'] = self.STATE_DEMOGRAPHICS

        # ----------------------------------------------------
        # 2. DEMOGRAPHICS INTAKE
        # ----------------------------------------------------
        elif current_state == self.STATE_DEMOGRAPHICS:
            age, gender = self.extract_demographics(user_message)
            if age is not None:
                session['age'] = age
            if gender is not None:
                session['gender'] = gender

            gender_label = "Erkek" if session['gender'] == 'M' else "Kadın"
            response_text = (
                f"Demografik profil kaydedildi: **{session['age']} yaşında, {gender_label}**.\n\n"
                "Özellikle akciğer ve solunum yolu riskini belirleyen en önemli faktörlerden biri tütün maruziyetidir:\n"
                "👉 **Sigara kullanıyor musunuz veya geçmişte kullandınız mı?**\n"
                "(Günde kaç paket ve kaç yıl boyunca? Örn: *'Günde 1 paket, 25 yıldır'* veya *'Hiç içmedim'*)"
            )
            quick_replies = [
                "Hiç içmedim (0 paket-yıl)",
                "Günde 1 paket, 10 yıl (10 paket-yıl)",
                "Günde 1 paket, 25 yıl (25 paket-yıl)",
                "Günde 2 paket, 30 yıl (60 paket-yıl)"
            ]
            session['state'] = self.STATE_SMOKING

        # ----------------------------------------------------
        # 3. SMOKING HISTORY
        # ----------------------------------------------------
        elif current_state == self.STATE_SMOKING:
            pack_years = self.extract_smoking_pack_years(user_message)
            session['smoking_pack_years'] = pack_years

            response_text = (
                f"Tütün maruziyeti hesaplandı: **{pack_years} Paket-Yıl**.\n\n"
                "Kalıtsal yatkınlık ve genetik mutasyon riskini değerlendirmek için:\n"
                "👉 **Ailenizde (1. derece akrabalarınız: anne, baba, kardeş veya çocuklarınızda) kanser öyküsü var mı?**"
            )
            quick_replies = ["Evet, ailemde kanser öyküsü var", "Hayır, ailemde kanser öyküsü yok"]
            session['state'] = self.STATE_FAMILY

        # ----------------------------------------------------
        # 4. FAMILY HISTORY
        # ----------------------------------------------------
        elif current_state == self.STATE_FAMILY:
            has_fam = self.extract_family_history(user_message)
            session['family_history'] = has_fam
            fam_str = "Var (1. Derece Pozitif)" if has_fam else "Yok (Negatif)"

            response_text = (
                f"Aile öyküsü kaydedildi: **{fam_str}**.\n\n"
                "Son olarak belirtilerinizin zamansal seyrini netleştirelim:\n"
                "👉 **Bu şikayetleriniz yaklaşık ne kadar süredir devam ediyor?**"
            )
            quick_replies = ["1 haftadan az", "2 - 3 hafta", "1 - 3 aydır", "6 aydan uzun süredir"]
            session['state'] = self.STATE_DURATION

        # ----------------------------------------------------
        # 5. DURATION & CLINICAL RISK EVALUATION
        # ----------------------------------------------------
        elif current_state == self.STATE_DURATION:
            session['duration_desc'] = user_message.strip()

            # Execute comprehensive calculation with CancerRiskPredictor
            eval_result = self.predictor.evaluate_patient(
                age=session['age'],
                gender=session['gender'],
                symptoms=session['symptoms'],
                smoking_pack_years=session['smoking_pack_years'],
                family_history=session['family_history']
            )

            p_info = eval_result['patient_profile']
            r_eval = eval_result['risk_evaluation']
            c_evid = eval_result['cohort_empirical_evidence']
            report_data = eval_result

            # Construct Comprehensive Clinical Consultation Report
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
                f"Teşekkür ederim. Verdiğiniz bilgiler **100.000 kişilik doğrulanmış klinik veri tabanımız** "
                "ve uluslararası kılavuzlar ile analiz edildi:\n\n"
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                f"📋 **KLİNİK ERKEN TEŞHİS DEĞERLENDİRME RAPORU**\n"
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                f"👤 **Hasta:** {p_info['age']} Yaşında {p_info['gender']} | Sigara: {p_info['smoking_pack_years']} paket-yıl | Aile Öyküsü: {'Var' if p_info['family_history'] else 'Yok'}\n"
                f"⏱️ **Süre:** {session['duration_desc']}\n"
                f"🚩 **Saptanan Belirtiler:** {', '.join(p_info['active_symptoms']) if p_info['active_symptoms'] else 'Akut Kırmızı Bayrak Yok'}\n\n"
                f"📊 **Bileşik Kanser Riski:** %{r_eval['composite_risk_percent']} -> **[{r_eval['risk_tier']}]**\n"
                f"{nice_alert}\n"
                f"🏥 **Öncelikli İncelenecek Alanlar:**\n{top_organs_text}\n"
                f"🔬 **Önerilen Öncelikli Test:** {r_eval['primary_diagnostic_recommendation']}\n\n"
                f"📈 **100.000 Kişilik Kohort Veri Eşleşmesi:**\n"
                f"  • Yaş Grubu ({c_evid['age_bracket']}): **{c_evid['matched_cohort_patients']:,}** benzer hasta incelendi.\n"
                f"  • Doğrulanmış Vaka: **{c_evid['cohort_diagnosed_cancers']:,}** kişide malignite saptandı (%{c_evid['cohort_incidence_percent']} insidans).\n\n"
                f"💡 **Hekim Önerisi & Eylem Planı:**\n{r_eval['clinical_action']}\n\n"
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
                "Yeni bir değerlendirme yapmak isterseniz *'Yeniden Başla'* yazabilirsiniz."
            )

            quick_replies = ["Yeniden Başla", "NICE Kılavuzu Nedir?", "Raporu İndir"]
            session['state'] = self.STATE_COMPLETED

        # ----------------------------------------------------
        # 6. COMPLETED OR FAQ STATE
        # ----------------------------------------------------
        elif current_state == self.STATE_COMPLETED:
            if re.search(r'\b(yeniden|baştan|tekrar|yeni|sıfırla)\b', user_message.lower()):
                return self.start_session(session_id)
            elif "nice" in user_message.lower():
                response_text = (
                    "ℹ️ **NICE NG12 Kılavuzu:** İngiltere Ulusal Sağlık ve Bakım Mükemmelliği Enstitüsü "
                    "(National Institute for Health and Care Excellence) tarafından yayımlanan kanser erken tanı rehberidir. "
                    "Herhangi bir belirtinin kanser pozitif prediktif değeri (PPV) %3 ve üzerine çıktığında, "
                    "hastanın 2 hafta içinde uzman tetkikine (2-Week-Wait) alınmasını zorunlu kılar."
                )
                quick_replies = ["Yeniden Başla", "Riskimi Tekrar Göster"]
            else:
                response_text = (
                    "Bu değerlendirme tamamlandı. Yeni bir şikayet veya hasta profili girmek için "
                    "aşağıdaki düğmeden **Yeniden Başla**'ya tıklayabilirsiniz."
                )
                quick_replies = ["Yeniden Başla"]

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
        print("=" * 70)
        print("🤖 ERKEN KANSER TEŞHİSİ & KLİNİK ANAMNEZ SORU BOTU (100k Kohort)")
        print("Çıkmak için 'q' veya 'exit' yazabilirsiniz.")
        print("=" * 70)

        init_data = self.start_session("cli_session")
        print(f"\n[BOT]: {init_data['message']}\n")
        if init_data.get('quick_replies'):
            print(f"  Hızlı Seçenekler: {', '.join(init_data['quick_replies'])}\n")

        while True:
            try:
                user_msg = input("Siz > ").strip()
                if not user_msg:
                    continue
                if user_msg.lower() in ['q', 'exit', 'quit', 'çıkış']:
                    print("\nSoru botu oturumu sonlandırıldı. Sağlıklı günler dileriz!\n")
                    break

                res = self.process_message(user_msg, "cli_session")
                print(f"\n[BOT]:\n{res['message']}\n")
                if res.get('quick_replies'):
                    print(f"  Hızlı Seçenekler: {', '.join(res['quick_replies'])}\n")

                if res.get('state') == self.STATE_COMPLETED:
                    print("-" * 70)

            except (KeyboardInterrupt, EOFError):
                print("\nSoru botu oturumu sonlandırıldı.\n")
                break


if __name__ == '__main__':
    bot = DiagnosticChatbot()
    bot.interactive_cli_session()
