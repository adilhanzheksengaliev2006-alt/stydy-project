
# stydy-project
# Биоинформатический анализ выживаемости и ландшафта мутаций при раке молочной железы (данные METABRIC)

**Автор:** Адильхан Жексенгалиев 
**Цель работы:** Исследовать частоту мутаций в когорте пациентов с раком молочной железы, оценить их совместную встретимость и выявить независимые генетические маркеры, влияющие на общую выживаемость.



## 1. Топ-10 самых часто мутирующих генов
На первом этапе через SQL-запрос к базе данных была посчитана частота встречаемости мутаций для каждого гена (исключая синонимичные замены `Silent`). Всего в анализируемую когорту вошел 1981 пациент.

<img width="1200" height="750" alt="top_mutated_genes" src="https://github.com/user-attachments/assets/0d17eb20-5c4a-4507-8613-4a4802719552" />

---

## 2. Попарный анализ взаимного исключения и синергии генов
Для выявления скрытых связей между путями канцерогенеза был проведен точный тест Фишера для всех возможных пар генов из Топ-10. Результаты визуализированы в виде тепловой карты Log Odds Ratio.

<img width="800" height="800" alt="mutation_co_occurrence_heatmap" src="https://github.com/user-attachments/assets/ff23a758-877d-4c1e-83e4-3909712e35a0" />


### Выводы по карте:
* **Взаимное исключение (Mutual Exclusivity):** Мутации в гене `TP53` практически не пересекаются с мутациями в `GATA3` ($p = 1.81 \times 10^{-24}$) и `PIK3CA` ($p = 2.56 \times 10^{-12}$). Это строгое математическое подтверждение существования двух независимых молекулярных путей развития опухоли (например, люминального и базальноподобного подтипов).
* **Совместная встречаемость (Co-occurrence):** Обнаружена значимая синергия между генами `PIK3CA` и `MAP3K1` (красный кластер), что указывает на их совместное участие и потенциальную кооперацию в онкогенезе.

---

## 3. Анализ выживаемости Каплана-Мейера
Для оценки изолированного влияния мутации ключевого онкосупрессора `TP53` на общую выживаемость был построен классический график Каплана-Мейера.

<img width="1200" height="900" alt="kaplan_meier_tp53" src="https://github.com/user-attachments/assets/b0063737-a4d3-4651-9ef3-90d66ada87e4" />



### Вывод:
Мутация в гене `TP53` статистически значимо ассоциирована с крайне неблагоприятным прогнозом ($p < 0.0001$). Кривая выживаемости пациентов с аберрациями в данном гене падает значительно быстрее, снижая медиану продолжительности жизни.

---

## 4. Многофакторный анализ (Модель регрессии Кокса)
Чтобы оценить вклад каждого гена с поправкой на присутствие всех остальных альтераций, была построена многофакторная модель пропорциональных рисков Кокса для всего Топ-10 генов одновременно.

<img width="800" height="700" alt="cox_forest_plot" src="https://github.com/user-attachments/assets/f0106655-b203-4257-aac3-b218007abaea" />



### Ключевые результаты модели:
* **Ген TP53 (HR = 1.24, p = 0.001):** Подтвержден как сильный независимый маркер неблагоприятного прогноза. Наличие мутации увеличивает риск смерти пациента на 24% вне зависимости от статуса других генов.
* **Ген GATA3 (HR = 0.62, p < 0.001):** Действует как прогностически благоприятный маркер, снижая риск смерти на 38%. Мутации в нем характерны для менее агрессивных форм рака.
* Остальные гены из Топ-10 (включая *PIK3CA* и *MUC16*) в рамках многофакторного анализа не показали статистически значимого независимого влияния на общую выживаемость ($p > 0.05$), выступая, вероятно, в роли генов-пассажиров.

* [Breast Cancer, METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) 
через cBioPortal — ~2000 пациентов с раком молочной железы, 
мутации + клинические данные о выживаемости


мои скрипт на языке R:

<details>
<summary>1. Подготовка данных</summary>

```r
library(RSQLite)

con <- dbConnect(SQLite(), "cancer_data.db")

clinical <- dbGetQuery(
  con,
  "SELECT PATIENT_ID, OS_MONTHS, OS_STATUS FROM clinical"
)

mutations <- dbGetQuery(
  con,
  "SELECT PATIENT_ID, Hugo_Symbol, Variant_Classification
   FROM mutations"
)

dbDisconnect(con)


top_genes <- c(
  "PIK3CA",
  "TP53",
  "MUC16",
  "AHNAK2",
  "SYNE1",
  "KMT2C",
  "GATA3",
  "MAP3K1",
  "CDH1",
  "DNAH11"
)


mutations <- mutations[
  mutations$Variant_Classification != "Silent",
]


data <- clinical


for (gene in top_genes) {

  patients <- unique(
    mutations$PATIENT_ID[
      mutations$Hugo_Symbol == gene
    ]
  )

  data[[gene]] <- ifelse(
    data$PATIENT_ID %in% patients,
    1,
    0
  )
}


data$TIME <- as.numeric(
  as.character(data$OS_MONTHS)
)


data$EVENT <- ifelse(
  grepl(
    "DECEASED|^1",
    as.character(data$OS_STATUS),
    ignore.case = TRUE
  ),
  1,
  0
)


data <- data[
  !is.na(data$TIME) &
  !is.na(data$EVENT),
]


head(data)






<summary>Развернуть R-код проекта</summary>

```r
ТУТ ВЕСЬ ТВОЙ R-КОД

ggforest(cox_model, data = data)
```
</details>

<details>
<summary>2. Анализ взаимодействия генов</summary>

```r
library(RSQLite)
library(corrplot)


con <- dbConnect(SQLite(), "cancer_data.db")

clinical <- dbGetQuery(
  con,
  "SELECT PATIENT_ID FROM clinical"
)

mutations <- dbGetQuery(
  con,
  "SELECT PATIENT_ID, Hugo_Symbol, Variant_Classification
   FROM mutations"
)

dbDisconnect(con)


top_genes <- c(
  "PIK3CA",
  "TP53",
  "MUC16",
  "AHNAK2",
  "SYNE1",
  "KMT2C",
  "GATA3",
  "MAP3K1",
  "CDH1",
  "DNAH11"
)


mutations <- mutations[
  mutations$Variant_Classification != "Silent",
]


data <- clinical


for (gene in top_genes) {

  patients <- unique(
    mutations$PATIENT_ID[
      mutations$Hugo_Symbol == gene
    ]
  )

  data[[gene]] <- ifelse(
    data$PATIENT_ID %in% patients,
    1,
    0
  )
}


n <- length(top_genes)

p_values <- matrix(
  1,
  nrow = n,
  ncol = n
)

odds_ratios <- matrix(
  0,
  nrow = n,
  ncol = n
)

rownames(p_values) <- top_genes
colnames(p_values) <- top_genes

rownames(odds_ratios) <- top_genes
colnames(odds_ratios) <- top_genes


for (i in 1:(n - 1)) {

  for (j in (i + 1):n) {

    gene1 <- data[[top_genes[i]]]
    gene2 <- data[[top_genes[j]]]

    table_genes <- table(
      factor(gene1, levels = c(0, 1)),
      factor(gene2, levels = c(0, 1))
    )

    test <- fisher.test(table_genes)

    odds_ratio <- test$estimate

    if (!is.finite(odds_ratio)) {
      odds_ratio <- 1
    }

    odds_ratios[i, j] <- log(odds_ratio)
    odds_ratios[j, i] <- log(odds_ratio)

    p_values[i, j] <- test$p.value
    p_values[j, i] <- test$p.value
  }
}


corrplot(
  odds_ratios,
  method = "color",
  is.corr = FALSE,
  p.mat = p_values,
  sig.level = 0.05,
  insig = "blank"
)
```
</details>

<details>
<summary>3. Анализ выживаемости</summary>

```r
library(RSQLite)
library(survival)
library(survminer)


con <- dbConnect(SQLite(), "cancer_data.db")

clinical <- dbGetQuery(
  con,
  "SELECT PATIENT_ID, OS_MONTHS, OS_STATUS
   FROM clinical"
)

mutations <- dbGetQuery(
  con,
  "SELECT PATIENT_ID, Hugo_Symbol, Variant_Classification
   FROM mutations"
)

dbDisconnect(con)


top_genes <- c(
  "PIK3CA",
  "TP53",
  "MUC16",
  "AHNAK2",
  "SYNE1",
  "KMT2C",
  "GATA3",
  "MAP3K1",
  "CDH1",
  "DNAH11"
)


mutations <- mutations[
  mutations$Variant_Classification != "Silent",
]


data <- clinical


for (gene in top_genes) {

  patients <- unique(
    mutations$PATIENT_ID[
      mutations$Hugo_Symbol == gene
    ]
  )

  data[[gene]] <- ifelse(
    data$PATIENT_ID %in% patients,
    1,
    0
  )
}


data$TIME <- as.numeric(
  as.character(data$OS_MONTHS)
)


data$EVENT <- ifelse(
  grepl(
    "DECEASED|^1",
    as.character(data$OS_STATUS),
    ignore.case = TRUE
  ),
  1,
  0
)


data <- data[
  !is.na(data$TIME) &
  !is.na(data$EVENT),
]


data$TP53_Status <- ifelse(
  data$TP53 == 1,
  "TP53 mutation",
  "No TP53 mutation"
)


fit <- survfit(
  Surv(TIME, EVENT) ~ TP53_Status,
  data = data
)


ggsurvplot(
  fit,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE
)


formula <- as.formula(
  paste(
    "Surv(TIME, EVENT) ~",
    paste(top_genes, collapse = " + ")
  )
)


cox_model <- coxph(
  formula,
  data = data
)


ggforest(
  cox_model,
  data = data
)
```
</details>


<details>
<summary>А вот мой SQL-код</summary>

```python
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt

MUTATIONS_FILE = "metabric_data/data_mutations.txt"
CLINICAL_FILE = "metabric_data/data_clinical_patient.txt"
GENE_OF_INTEREST = "TP53"

mutations = pd.read_csv(
    MUTATIONS_FILE,
    sep="\t",
    comment="#"
)

clinical = pd.read_csv(
    CLINICAL_FILE,
    sep="\t",
    comment="#"
)

patient_ids = clinical["PATIENT_ID"].tolist()

def find_patient(sample):
    for patient in patient_ids:
        if str(sample).startswith(patient):
            return patient
    return None

mutations["PATIENT_ID"] = mutations[
    "Tumor_Sample_Barcode"
].apply(find_patient)

conn = sqlite3.connect("cancer_data.db")

mutations.to_sql(
    "mutations",
    conn,
    if_exists="replace",
    index=False
)

clinical.to_sql(
    "clinical",
    conn,
    if_exists="replace",
    index=False
)

query_genes = """
SELECT
    Hugo_Symbol AS gene,
    COUNT(DISTINCT PATIENT_ID) AS patients
FROM mutations
GROUP BY Hugo_Symbol
ORDER BY patients DESC
LIMIT 15
"""

top_genes = pd.read_sql(
    query_genes,
    conn
)

print(top_genes)

query_survival = f"""
SELECT
    CASE
        WHEN mutations.PATIENT_ID IS NOT NULL
        THEN 'Mutation'
        ELSE 'No mutation'
    END AS group_name,

    COUNT(DISTINCT clinical.PATIENT_ID) AS patients,

    ROUND(
        AVG(clinical.OS_MONTHS),
        1
    ) AS average_survival

FROM clinical

LEFT JOIN (
    SELECT DISTINCT PATIENT_ID
    FROM mutations
    WHERE Hugo_Symbol = '{GENE_OF_INTEREST}'
) AS mutations

ON clinical.PATIENT_ID = mutations.PATIENT_ID

GROUP BY group_name
"""

survival = pd.read_sql(
    query_survival,
    conn
)

print(survival)

conn.close()

plt.figure(figsize=(8, 5))

plt.barh(
    top_genes["gene"],
    top_genes["patients"]
)

plt.xlabel("Patients")
plt.ylabel("Gene")
plt.title("Most frequently mutated genes")

plt.gca().invert_yaxis()

plt.tight_layout()

plt.savefig(
    "top_mutated_genes.png",
    dpi=150
)
```

</details>
