# Bilimsel Araştırma Raporu: Uçtan Uca Biyoinformatik ve İlaç Keşfi Hattı

Bu rapor, Antigravity **Science** araç seti kullanılarak gerçekleştirilen çok katmanlı moleküler biyoloji, genomik varyant analizi, farmakoloji ve bilimsel literatür taraması sonuçlarını içermektedir.

---

## 1. Giriş ve Biyolojik Çerçeve
Araştırma kapsamında insan eritrositlerinde oksijen taşınmasını yöneten **Hemoglobin** kompleksi ve dünyada en sık görülen monogenik hastalıklardan biri olan **Orak Hücreli Anemi (Sickle Cell Disease)** ele alınmıştır.

---

## 2. 3B Yapı ve Fonksiyonel Bölge Analizi
* **Protein Veritabanı (PDB):** `4HHB` (İnsan Deoksihemoglobini)
* **Alt Birimler:** Heterotetramer ($\alpha_1\beta_1\alpha_2\beta_2$)
* **İncelenen Kritik Kalıntılar:**
  * **Proksimal Histidin (`Chain A: His87`):** $Fe^{2+}$ demir atomuna doğrudan kovalent bağla tutunur.
  * **Distal Histidin (`Chain A: His58`):** Oksijen bağlanmasını stabilize eder, CO bağlanmasını engeller.
  * **Oraklaşma Mutasyon Noktası (`Chain B/D: Glu6`):** Mutasyon sonucu Valin kalıntısına dönüşerek deoksi formda lifsi polimerlerin oluşmasına ve alyuvar deformasyonuna yol açar.

---

## 3. Genomik Varyant ve Klinik Önem
* **Gen:** *HBB* (*Hemoglobin Subunit Beta*)
* **Kritik Varyant:** `rs334` (`chr11:5227002:T>A`, HGVS: `NM_000518.5:c.20A>T`, `p.Glu7Val` / `Glu6Val`)
* **ClinVar Sınıflandırması:** **Pathogenic (Patojenik)** (Variation ID: `15333`)
* **Kontrol Karşılaştırması:** `BRCA2` genindeki `rs28897700` (`c.198A>G`, `p.Gln66=`) eşanlamlı varyantı ClinVar tarafından **Benign (Zararsız)** olarak sınıflandırılmıştır (Uzman paneli onaylı).

---

## 4. Farmakoloji ve İlaç Keşfi
* **Hedef İlaç:** **Voxelotor (Oxbryta)** (PubChem CID: `71602803`)
* **Formül:** $C_{19}H_{19}N_3O_3$ | **Molekül Ağırlığı:** 337.4 g/mol
* **Etki Mekanizması:** Hemoglobinin $\alpha$-zincirinin N-terminal valin kalıntısına kovalent/geri-dönüşümlü bağlanır. Hemoglobinin $O_2$ afinitesini artırarak oksi-R durumunu stabilize eder ve HbS polimerleşmesini durdurur.
* **Metformin Analizi (ChEMBL: `CHEMBL1431`):** Tip 2 diyabet ve metabolik regülasyonda kullanılan küçük molekül standart profil çıkarımı yapılmıştır.

---

## 5. Hücresel Yolaklar (Reactome)
*HBB* geninin yer aldığı birincil biyokimyasal mekanizmalar:
1. `R-HSA-1247673`: *Erythrocytes take up oxygen and release carbon dioxide*
2. `R-HSA-1480926`: *O2/CO2 exchange in erythrocytes*
3. `R-HSA-9927020`: *Heme assimilation*

---

## 6. CRISPR-Cas9 Devrimi ve Son Literatür
* **Makale:** *"Exagamglogene Autotemcel for Severe Sickle Cell Disease"*
* **Dergi:** *The New England Journal of Medicine (NEJM)* (PMID: 38661449)
* **Klinik Çözüm:** CRISPR-Cas9 aracılığıyla *BCL11A* eritroid güçlendiricisinin susturulması, fetal hemoglobinin (**HbF**) yeniden üretimini tetikleyerek hastalığı fonksiyonel olarak tedavi eder.
