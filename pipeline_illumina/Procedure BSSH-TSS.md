# Procédure BSSH et TSS

Procédure pour le transfert des sorties de séquenceur Illumina vers les services BaseSpace Sequence Hub (BSSH) et TruSight (TSS).

0. Configuration
1. Définir les cas à partir du fichier de soumission
2. Transférer les résultats de la run sur BSSH:
    - Séquençage
    - Transfert des BCLs
    - Conversion/demux en FASTQs
3. Connecter les FASTQs à TSS
4. Créer les cas pour les analyses TSS
5. Archiver et nettoyer

No. | Étape                                              | Outil
====|====================================================|===================
1   | Définir les cas à partir du fichier de soumission  | `define_cases.py`
2.1 | Transfert des BCLs                                 | `SampleSheet.csv`
2.2 | Conversion/demux en FASTQs                         | `SampleSheet.csv`
3   | Connecter les FASTQs à TSS                         | `SampleSheet.csv`
4.1 | Créer les cas, définis à l'étape 1, dans TSS       | `define_cases.py`
4.2 | Analyser les cas dans TSS                          | 

## 0. Configuration

Préparer les pré-requis et configurer l'environnement de travail et les outils.

### Pré-requis

- `bs`: Outil en ligne de commande pour interagir avec BSSH. Suivre les instructions d'Illumina pour l'installation et la configuration. L'outil `bs` est déjà installé sur spxp-app02.
- Jetons d'authentification (tokens): Obtenir auprès d'Illumina, en suivant les instructions sur leur [site]()
- Fichiers de soumission: Définition des cas
- SampleSheet V2: Description des échantillons et paramètres de séquençage. 

### Fichier de configuration

Configuration file expected at "~/.illumina/gapp_conf.json" and JSON content should have the following format and information:

```json
{
    "instance"         : "cac1.trusight.illumina.com",
    "X-ILMN-Domain"    : "chusj",
    "X-ILMN-Workgroup" : "42948014-b206-320d-b304-1af26fc98af3",
    "X-Auth-Token"     : "APIKey *#u37t_5KmQ4FWGfBl)Y)1",
    "testDefinitionId" : "278b1d65-4cad-44e1-89d6-425c26564380",
    "phenotips_user"   : "USERNAME",
    "phenotips_pass"   : "PASSWORD"
}
```

## 1. Définir les cas

Les cas (Case) sont décrits dans le fichier de soumission. Ce dernier est se trouve sur spxp-app02 dans le répertoire `spxp-app02://mnt/vs_nas_chusj/PUBLIC_CRHSJ-centre_genomique/Laboratoire/Nanuq/Reception_echantillons`.

## 2. Transférer les résultats

### 2.1 Transférer des BCLs

La déconvolution-conversion des BCLs en FASTQs se fait automatiquement grâce aux instructions inscrites dans la SampleSheet V2. Pour information, la SampleSheet V2 peut être récupérée _via_ l'API Web de Nanuq avec la commande: `wget --post-data 'j_username=USERNAME&j_password=PASSWORD' --no-cookies "https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/sampleSheetV2/NovaSeq/RUNNAME/" -O "SampleSheet.csv"`. 

Normalement, les runs Rapidomics sont configurées de sorte à ce que les fichiers BCLs sont téléversés sur BSSH à la suite du séquençage. (Cependant, la SampleSheet doit être nommée `SampleSheet.csv`.) #TODO[2022-04-28]: Vérifier. Dans la situation où des FASTQ au lieu des BCLs doivent être téléversés sur BSSH, il faut utiliser l'utilitaire d'Illumina `tss-cli` (recommandé) comme décrit dans la section 2.1 ou l'ancien outil `wgs-cli` (_c.f._ section 2.2).

## 2.2 Transferts avec `tss-cli`

L'utilitaire CLI d'Illumina `tss-cli-<version>` permet de transférer les FASTQ en ligne de commande. Ce dernier peut être téléchargé depuis le lien _Download CLI Tool_ sous le menu _Account_ dans notre compte TSS (_c.f._ [Guide](https://support.illumina.com/help/TruSight_Software_Suite_1000000110864/Content/Source/Informatics/VariantSW/TSSS/CLI_swTSIS.htm#)). Copier le jar dans `${HOME}/bin/`.

Pré-requis:

1. Avant d'être utilisé, il faut configurer l'outil avec notre clé API tel que décrit sur [Aide en-ligne d'Ilumina pour importer des fastq avec le CLI pour TSS](https://support.illumina.com/help/TruSight_Software_Suite_1000000110864/Content/Source/Informatics/VariantSW/TSSS/ImportFiles_swTSIS.htm)

    ```bash
    java -jar ~/bin/tss-cli-<version>.jar configure --domain chusj --workgroup 12345-b172-321a-b456-1bc10ca11fr6 --url chusj.cac1.trusight.illumina.com --apiKey 'Y0uR_AP1%k3yG03sH3r3!'
    ```
2. Préparer un dossier local pour le staging du transfert. L'utilitaire `tss-cli` requiert que tous les fastq d'un échantillon soit contenu dans un même répertoire. En se basant sur les identifiants d'échantillons de la SampleSheet, regrouper tous les fastq dans un même répertoire pour le transfert.
3. L'outil `tss-cli` requiert un fichier "manifest" en format JSON, pour lui indiquer les détails des fichiers à transférer [[ref.](https://support-docs.illumina.com/SW/TSSS/TruSight_SW_Suite/Content/SW/TSSS/SampleManifestJSON.htm)]. Le script `tss/make_fastq_manifest.py` permet de créer un fichier JSON manifest à partir d'un dossier contenant les fastq.

```json
{
    "files": [
        {
            "lane": "1",
            "library": "TruSeq_DNA_PCR-free",
            "read1": "/scratch/hien/_fastq/14126/14126_S1_L001_R1_001.fastq.gz",
            "read2": "/scratch/hien/_fastq/14126/14126_S1_L001_R2_001.fastq.gz",
            "readGroup": "14126_L001"
        },
        {
            "lane": "2",
            "library": "TruSeq_DNA_PCR-free",
            "read1": "/scratch/hien/_fastq/14126/14126_S1_L002_R1_001.fastq.gz",
            "read2": "/scratch/hien/_fastq/14126/14126_S1_L002_R2_001.fastq.gz",
            "readGroup": "14126_L002"
        },
        {
            "lane": "3",
            "library": "TruSeq_DNA_PCR-free",
            "read1": "/scratch/hien/_fastq/14126/14126_S1_L003_R1_001.fastq.gz",
            "read2": "/scratch/hien/_fastq/14126/14126_S1_L003_R2_001.fastq.gz",
            "readGroup": "14126_L003"
        },
        {
            "lane": "4",
            "library": "TruSeq_DNA_PCR-free",
            "read1": "/scratch/hien/_fastq/14126/14126_S1_L004_R1_001.fastq.gz",
            "read2": "/scratch/hien/_fastq/14126/14126_S1_L004_R2_001.fastq.gz",
            "readGroup": "14126_L004"
        }
    ]
}
```

```bash
# cd ${SCRATCH}/_fastq/
module load java/11.0.2
for sample in 14048 14049 14050; do
    python ~/bin/tss/make_fastq_manifest.py ${SCRATCH}/_fastq/${sample} > ${SCRATCH}/_fastq/${sample}/${sample}.json
    java -jar ~/bin/tss-cli-2.1.0.jar samples create --manifest ${SCRATCH}/_fastq/${sample}/${sample}.json --verbose --sample-id ${sample} &
done
```

Il est possible d'ajouter l'option `--overwrite` si la commande retourne une erreur du type "sample already exists".

### 2.3a Transferts avec `wgs-cli` avec staging

**_N.B._** `tss-cli` est l'outil officiellement recommandé par _Illumina_ mais nous avons initialement rencontré des soucis avec. [`wgs-cli`](https://support.illumina.com/help/TruSight_Software_Suite_1000000110864/Content/Source/Informatics/VariantSW/TSSS/ImportFiles_swTSIS.htm) a été utilisé pour trnasférer quelques échantillons pendant le support technique d'_Illumina_ cherche à corriger les problèmes.

Afin de pouvoir réutiliser les même FASTQ pour plusieurs Biosamples, l'outil `wgs-cli` a besoin qu'on lui spécifie le nom d'un dossier "staging". 

```bash
# Use the stage command to create the directory and add files.
sample="14143"
java -jar ~/bin/wgs-cli-2.0.1.jar stage --stageDirectory staging_${sample} --localDirectory ${SCRATCH}/_fastq/${sample}
```

### 2.3b Transferts avec `wgs-cli` sans "staging"

Si les même FASTQs ne seront pas utilisés pour plusieurs échantillons, nous pouvons téléverser les FASTQ sans "staging" et les lier directement à un cas.

    When you upload a FASTQ file directly to a case, you directly link the file to a sample in the case. To upload FASTQ files that can be used for more than one case, [see Stage and Link FASTQ Files](https://support.illumina.com/help/TruSight_Software_Suite_1000000110864/Content/Source/Informatics/VariantSW/TSSS/ImportFiles_swTSIS.htm#Stage).

```bash
# Use the fastqs command to upload and associate the files with a sample
# java -jar wgs-cli.jar --configFilePath <configfilepath> fastqs --localDirectory <local directory path> --sampleId <sample id>
sample="14134" # mother=14135 father=14136
java -jar ~/bin/wgs-cli-2.0.1.jar fastqs --localDirectory ${SCRATCH}/_fastq/${sample} --sampleId ${sample}
```

**_N.B._** Si le CaseID avec les samples des trios n'est pas pralablement créé, la méthode `wgs-cli fastqs` rapporte une erreur du genre:

    2021-11-23T11:24:28.237045-05:00[EST5EDT]  Upload complete
    2021-11-23T11:24:28.902478-05:00[EST5EDT]  [ERROR]  Failure while attempting to associate fastqs with sample ID

Cependant les séquences sont bel et bien téléversées, et peuvent être retrouvées (cocher la case "link FASTQ" du Biosample) lors de la création du CaseID, si celle-ci est réalisée ultérieurement (_c.f._ étape "_Créer un CaseID_" ci-dessous).

## 3. Connecter les FASTQs à TSS

## 4. Créer les cas pour les analyses TSS

## 5. Archiver et nettoyer

Archiver les BCLs
Archiver les FASTQs (oui ou non?)
Vérifier l'intégrité des transferts
Supprimer les fichiers

## Exemples

### Créer un cas

Minimum de détail pour la création d'un MONO (Proband)

```json
{
  "comment": "Comment text",
  "displayId": "FOO-BAR-BAZ",
  "subjects": [
      {
          "gender": "FEMALE",
          "isAffected": "AFFECTED",
          "relationshipToProband": "PROBAND",
          "samples": [
              {
                  "externalSampleId": "12345",
                  "externalSampleName": "GMsample"
              }
          ],
          "reportTypes": [
              "05b3e363-298c-481e-a760-c2e457200069"
          ]
      }
  ],
  "tags": [
      "foo",
      "bar",
      "baz"
  ],
  "testDefinitionId": "278b1d65-4cad-44e1-89d6-425c26564380"
}
```

Et pour la création d'un DUO.

```json
{
  "comment": "Comment text",
  "displayId": "FOO-BAR-DUO",
  "subjects": [
    {
      "gender": "FEMALE",
      "isAffected": "AFFECTED",
      "relationshipToProband": "PROBAND",
      "samples": [{}],
      "reportTypes": ["05b3e363-298c-481e-a760-c2e457200069"]
    },
    {
      "gender": "MALE",
      "isAffected": "UNKNOWN",
      "relationshipToProband": "FATHER",
      "samples": [{}],
      "reportTypes": []
    }
  ],
  "tags": ["foo", "bar", "duo"],
  "testDefinitionId": "278b1d65-4cad-44e1-89d6-425c26564380"
}
```

## Références

- [TSS Web API](https://chusj.cac1.trusight.illumina.com/crs/swagger-ui/index.html#/Cases/createCaseUsingPOST)
- Emplacement des fichier de soumission: spxp-app02://mnt/vs_nas_chusj/PUBLIC_CRHSJ-centre_genomique/Laboratoire/Nanuq/Reception_echantillons

## Annexe

Extrait du courriel d'illumina qui résume les grandes lignes de cette procédure.

    From: Cutler, Kyle <kcutler@illumina.com>
    Sent: March 16, 2022 16:19
    Subject: RE: Trusight Software 2.6 Demo

    (...)
    1. API Links format https://<domain>.<region>.trusight.illumina.com/<service>/swagger-ui/index.html
        a. Example: Case Registry Service, https://ilmn-cci.cac1.trusight.illumina.com/crs/swagger-ui/index.html
        b. Test Definitions Service, TMS
        c. Filter Management Service, FMS
        d. Variant Query Service, VQS
        e. Knowledge Network Service, KNS
        f. IGV Data Visualization Service, DVS
        g. Draft Report Service, DRS
    2.Gene/Transcript List differences across genome build and transcript source
        a. Currently working with the team to share that full list of gene names, should be next week hopefully
    3. Automated Steps for Flow from Sequencer > BSSH > TSS
        a. Required components BSSH CLI, https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
        b. Upload Run to BSSH with V2 sample sheet required cloud data and cloud settings (Specifying sample IDs, ProjectName, and BSSH App to run BCL Convert)

            [Cloud_Settings]              
            BsshApp,Illumina-inc.bcl-convert.1.3.0
            [Cloud_Data]                                                                                               
            Sample_ID,DataAggregationGroup,ProjectName,LibraryName,LibraryPrepKitUrn,LibraryPrepKitName,IndexAdapterKitUrn,IndexAdapterKitName
            Test1_NA12878_IDPF_A01,,,,,,,
            Test1_IDPF_S4_24plex_110122,,,,,,,
        
        c. Use biosample manifest file to create biosamples in BSSH that will automatically launch TSS Connect (yield trigger optional)
            i. Using BSSH CLI example for a single sample:
            ii. bs -c <bssh_profile> biosample create -p "<project_name>" -n "<sample_id>" --analysis-workflow "TSS Connect v1.1" --prep-request="Unknown" --required-yield=90 --allow-existing
            iii. Upload biosample manifest file using “TSS Connect v1.1” in the analysis workflow field
            iv. https://support.illumina.com/help/BaseSpace_Sequence_Hub/Source/Informatics/BS/BiosampleWorkflowFiles_swBS.htm
        d. Create case using API or UI
            i. API using CRS https://ilmn-cci.cac1.trusight.illumina.com/crs/swagger-ui/index.html and POST /crs/api/v1/cases/ and specifying “externalSampleID” as specifying in previous SampleID steps, the only required fields are Test ID (from test definition) and associated report type (from test definition)
            ii. Create using UI, if sample ID entered in Case Creation UI after data is pushed, FASTQs will automatically be associated with the sample
        e. Optionally the TSS Connect BSSH App can be run independently of the automated method to push FASTQs/Biosamples to TSS
