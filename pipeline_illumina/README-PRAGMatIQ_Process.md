# Procédure PRAGatIQ

Cette procédure décrit comment créer des cas (trio, solo, duo ou quad) à analyser sur la plateforme d'[Emedgene (EMG)](https://chusaintejustine.emedgene.com/) une fois que les séquences des échantillons sont disponibles sur BaseSpace Sequence Hub (BSSH). La progression des opérations de séquençage puis de la déconvolution-conversion de la _Run_ peuvent être suivid sur BSSH sous l'onglet ["Run"](https://chusj.cac1.sh.basespace.illumina.com/runs/active) et ["Analyses"](https://chusj.cac1.sh.basespace.illumina.com/analyses/), respectivement.

En résumé, voici les étapes à suivre:

0. Préparer la création des analyses
    1. Obtenir les identifiants de connexion
    2. Créer un journal pour le suivi, nommé `README-${FC_SHORT}.ipynb`
    3. Se connecter à spxp-app02 `ssh ${USER}@10.128.80.26`
    4. Mettre en place l'environnement de travail `conda activate CQGC-utils`.
1. Récupérer les informations sur les familles dans Nanuq `python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_make_batch_from_nanuq.py ${FC_SHORT}`. Ce script génère le fichier CSV d'entrée emg_batch_manifest.csv, avec les chemins d'accès aux FASTQs sur BSSH.
2. Créer les cas sur Emedgene et lancer les analyses grâce au fichier CSV.
3. (**TODO**) Ajouter les participants _via_ l'API
4. Archiver les résultats

Où ${FC_SHORT} est le nom court de la _FlowCell/Run_. Ex: Si la _flowcell/Run_ se nomme "230727_A00516_0441_AHKVFYDMXY", ${FC_SHORT} est "A00516_0441".


## Procédure


### 0. Préparer la création des analyses

#### 1. Obtenir les identifiants de connexion

Demander à l'administrateur de fournir les droits d'accès (identifiant, mot-de-passe, permissions) aux services suivants:

- Nanuq
- Emedgene
- Phenotips
- BSSH

Pour Nanuq, il faut créer le fichier `~/.nanuq`, contenant une seule et unique ligne comme ceci:

    j_username=USER&j_password=PASS&toto=1

Remplacer 'USER' et 'PASS' par vos identifiants. Par exemple, avec la commande: `echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`

Pour les services Emedgene, Phenotips et BSSH, il faut créer un fichier texte de configuration nommé `~/.illumina/gapp_conf.json` dont le contenu ressemble à ceci:

    {
        "instance"         : "cac1.trusight.illumina.com",
        "X-ILMN-Domain"    : "chusj",
        "X-ILMN-Workgroup" : "TOKEN",
        "X-Auth-Token"     : "APIKey <TOKEN>",
        "testDefinitionId" : "TOKEN",
        "bs_apiServer"     : "https://api.cac1.sh.basespace.illumina.com",
        "bs_accessToken"   : "TOKEN",
        "X-Gene42-Server"  : "https://chusj.phenotips.com",
        "X-Gene42-Auth"    : "Basic <TOKEN>",
        "X-Gene42-Secret"  : "TOKEN"
    }

Un exemple de fichier est disponible sur [l'espace PrivateDoc dans GitHub](https://github.com/CQGC-Ste-Justine/PrivateDoc/blob/0f59d674fd5d91ca5da7880ab55cd753ad324203/gapp_conf.json).


#### 2. Créer un journal pour le suivi

Créer un journal pour le suivi, nommé `README-${FC_SHORT}.ipynb` et y inscrire en titre les noms de la _flowcell_ (ex.: "230711_A00516_0433_BHKVMJDMXY") et de l'exérience ("Seq_S2_PRAG_20230711"). Ces informations sont normalement communiqués dans un courriel du laboratoire CQGC, mais peuvent aussi être récupérées depuis BSSH (sous l'onglet "_Runs_").

Le format `ipynb` (Jupyter Notebook) est utilisé car du code peut y être exécuté, mais un simple format `txt` ou `md` peut aussi bien servir.

#### 3. Se connecter à spxp-app02 

`ssh ${USER}@10.128.80.26`

#### 4. Mettre en place l'environnement de travail 

`conda activate DEV-bcl_convert`.

TODO: `conda activate CQGC-utils && conda update pandas`


### 1. Récupérer les informations sur les familles

Les informations sur la constitution des familles ("Case") sont centralisées dans Nanuq. La commande suivante permet de récupérer la liste des échantillons sur la _Run_ et de fournir les données nécessaires à la création de cas dans Emedgene, incluant les termes HPOs associées aux patients dans la base de données [Phenotips](https://chusj-phenotips.us.auth0.com/). 

`python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_make_batch_from_nanuq.py --run ${FC_SHORT}`

La commande ci-dessus devrait générer une sortie d'écran comme ceci:

[2023-07-13@17:39:50] Cases for A00516_0433:

    sample_name biosample relation  gender labels      mrn cohort_type  \
    3  23-03938-T1     21939  PROBAND  FEMALE   CHUS  1538614        TRIO
    4  23-05449-T1     21940      MTH  FEMALE   CHUS  1279938        TRIO
    5  23-05448-T1     21941      FTH    MALE   CHUS   147112        TRIO
    0  23-05328-T1     21933  PROBAND    MALE   CHUS  1636084        TRIO
    1  23-05367-T1     21934      MTH  FEMALE   CHUS  1635873        TRIO
    2  23-05368-T1     21935      FTH    MALE   CHUS  1636411        TRIO

    date_of_birth(YYYY-MM-DD) phenotypes    family case_group_number  \
    3                12/10/2018        AFF  23-03938          P0000137
    4                02/10/1979        UNK  23-03938                na
    5                01/11/1977        UNK  23-03938                na
    0                16/06/2023        AFF  23-05328          P0000134
    1                25/12/1991        UNK  23-05328                na
    2                25/01/1990        UNK  23-05328                na

                                                                            hpos
    3             HP:0001251,HP:0001263,HP:0002064,HP:0011342,HP:0012447,HP:0030891
    4                                                                            na
    5                                                                            na
    0  HP:0000126,HP:0000201,HP:0000294,HP:0000316,HP:0000586,HP:0001631,HP:0012368
    1                                                                            na
    2                                                                            na

**_N.b._**: Les termes HPOs [Human Phenotype Ontology](https://hpo.jax.org/app/) sont contenus dans la base de données patients, Phenotips, et peuvent être consultées _via_ l'API avec un identifiant externe composé du "site + MRN" (ex. "CHUS1636084"). C'est ce que fait le script `emg_make_batch_from_nanuq.py` car l'identifiant primaire PhenotipsID (PID), n'est pas contenu dans Nanuq.

Alternativement, les informations sur les familles pour la _Run_, incluant le PID, peuvent être récupérées à partir du fichier Excel partagé sur "\\shsjcifs01\public crhsj\centre genomique\Projets\Pragmatiq\Étude PRAGMatIQ - Base de données.xlsx" ou sur Teams/OneDrive "Étude PRAGMatIQ - Base de données.xlsx" (canal de l'équipe _CHUSJ - Étude PRAGMatIQ_ dans Teams ("Documents > Général")). Ceci peut s'avérer nécessaire dans l'éventualité où le script `emg_make_batch_from_nanuq.py` échoue. _C.f._ Annexe A (Rechercher les termes HPOs) plus bas. En dernier recours, contacter Camille VARIN-TREMBLAY pour obtenir les informations manquantes.


### 2. Créer les cas et lancer les analyses

Le script Python `emg_make_batch_from_nanuq.py` génère automatiquement le fichier manifeste (CSV) d'entrée "emg_batch_manifest.csv", avec les chemins d'accès complets aux FASTQs sur BSSH. 

`python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_make_batch_from_nanuq.py ${FC_SHORT}`

Pour plus d'information au sujet du fichier d'entrée pour la création des cas par lot, consulter les spécifications d'Emedgene [CSV manifest specification](https://help.emedgene.com/en/articles/7231644-csv-format-requirements).

Puis se connecter _via_ l'interface web à Emedgene pour lancer les analyses.

1. Se connecter à [Emedgene](https://chusaintejustine.emedgene.com/) avec le compte cqgc.bioinfo.hsj@ssss.gouv.qc.ca.
2. Cliquer sur "+ New Case", en haut de la page à droite
3. Cliquer sur "Switch to batch"
4. Glisser-déposer le fichier CSV manifest dans la boîte "Browse files from your computer or drag it here ('csv' | Maximum file size: 5 MB, 50 cases)"
5. Cliquer sur "Launch Case".

Pour plus d'informations, voir les instructions d'_Emedgene_ sur comment [créer des cas par lot](https://r4a56nl8uxkx3w3a292kjabk.emedgene.com/articles/7221986-batch-case-upload). https://chusaintejustine.emedgene.com/v2/#/help/


### 4. (**TODO**) Ajouter les participants


### 5. Archiver les résultats

Une fois que les analyses sont terminées, rapattrier et archiver les résultats.

```bash
#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
cd ${HOME}/projects/ctb-rallard/COMMUN/PRAGMatIQ-EMG
#for sample in ${@}; do
for sample in GM231339 GM231362 GM231511 GM231395 GM231398 GM231399 GM231401 GM231455 GM231456 GM231429 GM231433 GM231434 GM231438 GM231437 GM231436 GM231435 GM231463 GM231471 GM231500 23-05696-T1 23-05748-T1 23-05746-T1 GM231495 GM231504 GM231505; do
    aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/${sample}/ ./${sample} --recursive
    echo ${sample} $( date +'%Y-%m-%d %T' ) >> archive_PRAGMatIQ.log
done
```


## Références

- [Production environment](https://chusaintejustine.emedgene.com)
- [Eval environment](https://stejustine.emedgene.com/)
- [Phenotips](https://chusj-phenotips.us.auth0.com/)
- [Phenotips API](https://docs.phenotips.com/reference/accessandauthentication)
- [Human Phenotype Ontology (HPOs)](https://hpo.jax.org/app/)
- [Batch case upload instructions by Emedgene](https://r4a56nl8uxkx3w3a292kjabk.emedgene.com/articles/7221986-batch-case-upload) 
- https://chusaintejustine.emedgene.com/v2/#/help/
- [CSV manifest specification](https://help.emedgene.com/en/articles/7231644-csv-format-requirements)
- Old/Alternative method, using scripts:
    - "Case_creation_script_v2.docx"
    - `create_batch_cases_v2.py --help`
    - "create_batch_cases_v2_template.csv"

## Annexe

### A. Rechercher les termes HPOs

Dans la situation où les PID doivent être récupérés manuellement pour chaque cas, utiliser le script `/staging2/soft/CQGC-utils/Analysis.pipeline_illumina/get_phenotips_by_id.py`. Par exemple:

```bash
python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/get_phenotips_by_id.py P0000134
# {'id': 'HP:0000126', 'label': 'Hydronephrosis'}
# {'id': 'HP:0000201', 'label': 'Pierre-Robin sequence'}
# {'id': 'HP:0000294', 'label': 'Low anterior hairline'}
# {'id': 'HP:0000316', 'label': 'Hypertelorism'}
# {'id': 'HP:0000586', 'label': 'Shallow orbits'}
# {'id': 'HP:0001631', 'label': 'Atrial septal defect'}
# {'id': 'HP:0012368', 'label': 'Flat face'}
#
#  HP:0000126, HP:0000201, HP:0000294, HP:0000316, HP:0000586, HP:0001631, HP:0012368
#
#  Hydronephrosis, Pierre-Robin sequence, Low anterior hairline, Hypertelorism, Shallow orbits, Atrial septal defect, Flat face
```

Ou directement depuis un Notebook `ipynb`, copier-coller le code suivant dans une cellule

```python
import sys
# sys.path.append("D:\\HSJ\\Workspace\\tss")
sys.path.append("/staging2/soft/CQGC-utils/lib/")
from gapp import Phenotips

pho = Phenotips()

def get_phenotips_hpos(phenotips_id):
    pids = pho.get_hpo(phenotips_id)
    ids    = []
    labels = []
    for pid in pids:
        labels.append(pid['label'])
        ids.append(pid['id'])
        print(pid)
    print("\n", ', '.join(ids))
    print("\n", ', '.join(labels))

# Replace the PID in the function's argument
get_phenotips_hpos(P0000134)
```

Créer les cas avec la commande (ne fonctionne pas bien). Privilégier l'interface web "Switch to batch" pour y déposer le fichier manifeste. `python create_batch_cases_v2.py -i emg_batch_manifest.csv -s 10123 -hu stejustine.emedgene.com -u cqgc.bioinfo.hsj@ssss.gouv.qc.ca -p 7TmbuM3TUCMwP -b`.


### B. Créer les cas par l'interface web

Pour créer les cas en utilisant l'interface web, il faut au préalable récupérer les informations sur les familles contenues dans le fichier Excel "Étude PRAGMatIQ - Base de données.xlsx" (canal de l'équipe _CHUSJ - Étude PRAGMatIQ_ dans Teams ("Documents > Général")) ou obtenues par le script `emg_make_batch_from_nanuq.py`.

1. Se connecter à [Emedgene](https://chusaintejustine.emedgene.com) avec l'identifiant 'cqgc.bioinfo.hsj@ssss.gouv.qc.ca'.
2. Cliquer sur le bouton "+ New Case", en haut à droite
    1. Dans le panneau de gauche, "Select sample type", choisir le bouton-radio "FASTQ/BAM/CRAM". Nous utiliserons les FASTQ, qui sont déconvolués automatiquement après la run de séquençage sur BaseSpace.
    2. Cliquer sur "Next", en bas à droite. Une nouvelle page nous permet de "Create family tree" (panneau de gauche) en renseignant chaque membre de la famille avec "Add patient info" (panneau de droite).
        1. Dans la boîte de texte "Clincal Notes" du panneau de gauche "Create family tree", copier-coller le PhenotipsID
3. Renseigner les informations du Proband
    1. "Add sample" > sélectionner le bouton-radio "Create new sample".
    2. Cliquer sur "Continue"
    3. "Upload new sample", 
        1. Copier-coller le nom de l'échantillon dans la boîte de texte "Sample Name:" 
        2. Sélectionner les FASTQs: 
            1. Cliquer sur le bouton "Choose from storage", Une fenêtre flottante apparaît pour sélectioner les ficiers
            2. Dans la liste déroulante de la source ("S3 - par défaut"), dérouler et choisir "BASESPACE - 69f74cfd101e4befaf8d54f04b057e27"
            3. Dans la boîte de sélection en dessous de la liste déroulante, naviguer le chemin d'accès jusqu'à `projects/PRAGMatIQ_{site}/biosamples/{biosample_id}/datasets/{biosample_id}_L{1-4}/sequenced files/` (remplacer {site}, {biosample_id}, et lanes L{1-4} par la valeur correspondante)
            4. Cliquer sur tous les fichiers `*.fastq.gz` présents dans le dossier pour les sélectionner
            5. Cliquer sur "Continue"
            6. Répéter les étapes de la section 3.3.2 pour tous les FASTQs, dans toutes les lanes (2 ou 4, selon le type de flowcell utilisé pour le séquençage)
    4. "Patient info" (étape 6, d'après l'interface d'Emedgene). Basé sur les informations récupérées de Nanuq ou de la base de données Excel:
        1. Choisr le "Gender" (Male|Female|Unknown)
        2. Renseigner la date de naissance mm/dd/yyyy (**_n.b._**: les dates dans nos tableaux sont en format dd/mm/yyy)
        3. "Proband Phenotypes", cliquer sur "Batch mode" pour pouvoir copier-coller plusieurs termes HPO à la fois
            1. Copier-coller la liste des termes HPO dans cette boîte de texte
        4. Cliquer sur "Complete", en vérifiant les informations. Prendre soin que les noms des FASTQs correspondent tous au même biosample, incluant tous les lanes et reads forward et reverse. Au besoin, cliquer sur "Edit" et corriger les erreurs.
4. Renseigner les informations des autres membres de la famille
    1. Dans le panneau de gauche , "Create family tree", cliquer sur "Proband"
    2. Cliquer sur le petit cercle "+"
    3. Dans le panneau à droite "Add patient information", répéter les étapes 3.1 à 3.3.6 ci-dessus.
        1. Cliquer sur "Complete", en vérifiant les informations. **_N.b._**: prendre soin que les noms des FASTQs correspondent tous au même biosample, incluant tous les lanes et reads forward et reverse. Au besoin, cliquer sur "Edit" et corriger les erreurs.
    4. Répéter les étapes 4.1 à 4.3 pour chaque membre de la famille. **_N.b._**: Bien cliquer sur l'icône du "Proband" à chaque fois qu'il faut créer un parent pour celui-ci. Sinon, il y a le risque d'ajouter un parent au parent.
5. Cliquer sur "Next", Une nouvelle page s'affiche avec un panneau "Case info" à gauche. 
6. Renseigner les détails de l'analyse dans le nouveau panneau de gauche "Case info":
    1. "Select case type", dans le menu déroulant "Case type", choisir "Whole Genome"
    2. "Sequencing information", choisir le bouton-radio "No kit" puis
        1. Cliquer sur "Continue"
    4. (Oui, 4, Emedgene ne sait pas compter :-P) "Select genes list", choisir le bouton-radio "All genes"
    7. "Select preset": laisser sur "Default"
    8. O"ptional: Additional case info"
        1. Copier-coller le PhenotipsID dans la boite de texte "Indication for testing"
        2. "Label", dans le menu déroulant, choisir le label correspondant au site qui a demandé l'analyse (CHUSJ, CHUS,...)
        3. Cliquer sur "Complete"
9. Ajouter les participants à la liste (_c.f._ cid-essous)
    1. Cliquer sur "Send".

    

### aws

Some alignement and quality metric files can be downloaded directly from the web interface, under the "Lab" tab.
Additional logs and output files (CRAM, VCF) can be accessed from Emedgene's AWS S3 bucketunder 'emg-auto-results/<customer>/'. For example:

```bash
# Assuming that "emedgene" profile for connecetion is already set up
cd ~/tmp
aws s3 --profile emedgene ls --recursive --page-size 100000 s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/GM230424
aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/GM230424/GM230424_vlocal_2023-02-28-07-04_sample.log .
grep "Uniformity" GM230424_vlocal_2023-02-28-07-04_sample.log
grep "uniformity" GM230424_vlocal_2023-02-28-07-04_sample.log
```

[Production environment](https://chusaintejustine.emedgene.com)

    [11:15] Silver, Talia
    You won't need boost genes, gene list ID, kit ID for your genomes

    [11:24] Silver, Talia
    Label IDs:

    0: {id: 3, name: "Type of cases", used: false}
    1: {id: 4, name: "NICU", used: false}
    2: {id: 5, name: "PICU", used: false}
    3: {id: 6, name: "Paediatric Wards", used: false}
    4: {id: 7, name: "Other", used: false}
    5: {id: 8, name: "Re-analysis driven by clinician", used: false}
    6: {id: 9, name: "Re-analysis 1 year", used: false}
    7: {id: 10, name: "Validation", used: true}
    8: {id: 11, name: "Do not charge", used: false}
    9: {id: 12, name: "CHUS", used: true}
    10: {id: 13, name: "CHUSJ", used: true}
    11: {id: 14, name: "CHUQ", used: false}
    12: {id: 15, name: "MUHC", used: false}

    Use the IDs inside the brackets

    [11:25] Silver, Talia
    Storage_id: 10126
    name: "Illumina basespace"
    CHU_Sainte_Justine

[Eval environment](stejustine.emedgene.com):

    [11:27] Silver, Talia
    Storage_id: 10123
    name: "Illumina basespace"
    Ste_Justine_Eval

    [11:47] 
    11:47 La réunion est terminée : 47 min 12 s


    From: Silver, Talia <tsilver@illumina.com>
    Sent: December 13, 2022 16:19
    To: Quang-Hien Le (HSJ) <quang-hien.le.hsj@ssss.gouv.qc.ca>; Daoud, Hussein <hdaoud@illumina.com>; Alexandre Dionne-Laporte (HSJ) <alexandre.dionne-laporte.hsj@ssss.gouv.qc.ca>; Rene Allard (HSJ) <rene.allard.hsj@ssss.gouv.qc.ca>; Fischer, Yair <yfischer@illumina.com>; Bouslama, Sidki <sbouslama@illumina.com>
    Cc: Faucher, David <dfaucher@illumina.com>; Camille Varin-Tremblay (HSJ) <camille.varin-tremblay.hsj@ssss.gouv.qc.ca>
    Subject: Re: CHUSJ EMG/TSS

    Hi Hien,

    First off, make sure the user that is running the script is logged in and is moved to Ste_Justine_eval.
    Then use the url: stejustine.emedgene.com
    This BSSH Storage ID: 10123

    These are the label ids for Ste_Justine_eval:

    {id: 1, name: "Validation", used: true}
    1: {id: 2, name: "Do Not Charge", used: true}
    2: {id: 3, name: "Case Re-Analysis", used: false}
    3: {id: 6, name: "Type of cases", used: false}
    4: {id: 7, name: "NICU", used: false}
    5: {id: 8, name: "PICU", used: false}
    6: {id: 9, name: "Paediatric Wards", used: false}
    7: {id: 10, name: "Other", used: false}
    8: {id: 11, name: "Re-analysis driven by clinician", used: false}
    9: {id: 12, name: "Re-analysis 1 year", used: false}
    10: {id: 13, name: "Do not charge", used: false}
    11: {id: 14, name: "CHUS", used: false}
    12: {id: 15, name: "CHUSJ", used: false}
    13: {id: 16, name: "CHUQ", used: false}
    14: {id: 17, name: "MUHC", used: false}
