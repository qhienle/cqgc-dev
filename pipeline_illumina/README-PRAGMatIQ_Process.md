# Procédure PRAGMatIQ

Cette procédure décrit comment créer des cas (trio, solo, duo ou quad) à analyser sur la plateforme d'[Emedgene (EMG)](https://chusaintejustine.emedgene.com/) une fois que les séquences des échantillons sont disponibles sur BaseSpace Sequence Hub (BSSH). Les BCLs issus des séquenceurs sont d'abord convertis en FASTQs par DRAGEN BCLConvert, soit directement sur le NovaSeqX, soit par notre serveur local ou sur BSSH. La progression des opérations de séquençage puis de la déconvolution-conversion de la _Run_ peuvent être suivies sur BSSH sous l'onglet ["Run"](https://chusj.cac1.sh.basespace.illumina.com/runs/active) et ["Analyses"](https://chusj.cac1.sh.basespace.illumina.com/analyses/), respectivement. Les cas à analyser sont finalement créés dans Emedgene à l'aide d'un fichier "batch manifest" en format CSV. Ce fichier contient les informations sur les individus des familles (sexe, âge, relation, FASTQs associés, termes HPOs, _etc_.). 

En résumé, voici les étapes à suivre:

0. Pré-requis
    1. Obtenir les identifiants de connexion
    2. Créer un journal pour le suivi, par exemple `README-${FC_SHORT}.ipynb` (*)
    3. Se connecter à spxp-app02 `ssh ${USER}@10.128.80.26`
    4. Mettre en place l'environnement de travail `conda activate CQGC-utils`.
1. Collecter les informations sur les familles et les métriques DragenGermline `emg_collect_dragen_metrics.py ${FC}`.
2. Téléverser les FASTQs sur BaseSpace `emg_upload_fastqs.py`
3. Créer les cas sur Emedgene 
    1. Générer le fichier "emg_batch_manifest.csv" `emg_make_batch_from_nanuq.py ${FC}`
    2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
    3. (**TODO**) Ajouter les participants _via_ l'API
4. Archiver les résultats
5. Nettoyer

Où ${FC_SHORT} est le nom court de la _FlowCell/Run_. Ex: Si la _flowcell/Run_ se nomme "230727_A00516_0441_AHKVFYDMXY", ${FC_SHORT} est "A00516_0441".

(*) Un fichier `README-FC_SHORT-template.ipynb` est disponible sur GitHub. Il contient le code à copier-coller et les cases à remplir pour le suivi des opérations.

***_N.B._***: Afin de connecter les informations génétiques des familles (issues de Nanuq) aux phénotypes (termes HPOs contenus dans Phenotips), il est impératif que les deux champs EP (Établissement Public) et MRN (Medical Record Number) soient bien renseignés dans les deux systèmes par les collègues. Sinon, il faut obtenir les informations de Phenotips (identifiants PID) par "Camille Varin-Tremblay (HSJ)" <camille.varin-tremblay.hsj@ssss.gouv.qc.ca>.


## Procédure


### 0. Pré-requis

Installation des outils et obtention des identifiants de connexion aux différents services.

#### 1. Obtenir les identifiants de connexion

Demander à l'administrateur de fournir les droits d'accès (identifiant, mot-de-passe, permissions) aux services suivants:

- ENG: Emedgene, pour les analyses des cas
- Nanuq: Base de données des échantillons de séquençage du CQGC
- Phenotips: Base de données contenant les termes HPO des patients
- BSSH: _i.e._ BaseSpace, stockage infonuagique d'Illumina où les FASTQs seront mis à disposition pour Emedgene
    - `bs`: Utilitaire CLI pour interagir avec BaseSpace, notamment pour téléverser les FASTQs depuis spxp-app02

Pour Nanuq, il faut créer le fichier `~/.nanuq`, contenant une seule et unique ligne comme ceci:

    j_username=USER&j_password=PASS&toto=1

Remplacer 'USER' et 'PASS' par vos identifiants. Par exemple, avec la commande: `echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`

Pour les services Emedgene, Phenotips et BSSH, il faut créer un fichier texte de configuration nommé `~/.illumina/gapp_conf.json` dont le contenu ressemble à ceci:

    {
        "instance"         : "cac1.trusight.illumina.com",
        "X-ILMN-Domain"    : "chusj",
        "X-ILMN-Workgroup" : "<replace_with_your_hash_key>",
        "X-Auth-Token"     : "APIKey <replace_with_your_hash_key>",
        "testDefinitionId" : "<replace_with_your_hash_key>",
        "bs_apiServer"     : "https://api.cac1.sh.basespace.illumina.com",
        "bs_accessToken"   : "<replace_with_your_hash_key>",
        "X-Gene42-Server"  : "https://chusj.phenotips.com",
        "X-Gene42-Auth"    : "Basic replace_with_your_hash_key",
        "X-Gene42-Secret"  : "<replace_with_your_hash_key>",
        "REDCap-Server"    : "https://tacc-redcap.bic.mni.mcgill.ca/api/",
        "REDCap-Token"     : "<replace_with_your_hash_key>"
    }

Un exemple de fichier est disponible sur [l'espace PrivateDoc dans GitHub](https://github.com/CQGC-Ste-Justine/PrivateDoc/blob/0f59d674fd5d91ca5da7880ab55cd753ad324203/gapp_conf.json).

De plus, l'utilitaire CLI `bs`, installé sur spxp-app02, va rechercher un fichier de configuration `*.cfg` dans votre "home", `/home/hienle/.basespace`. Voici un exemple de `spxp-app02://home/hienle/.basespace/cac1.cfg`

    apiServer   = https://api.cac1.sh.basespace.illumina.com
    accessToken = <TOKEN>

Veuillez vous référer à la documentation d'Illumina afin de générer le jeton d'accès (token).

#### 2. Créer un journal pour le suivi

Créer un journal dans lequel seront notés le suivi des opérations. Par exemple, dans un notebook Jupyter nommé `README-${FC_SHORT}.ipynb` et y inscrire en titre les noms de la _flowcell_ (_e.g._.: "230711_A00516_0433_BHKVMJDMXY", où `${FC_SHORT}` serait dans ce cas "A00516_0433") et de l'exérience ("Seq_S2_PRAG_20230711"). Ces informations sont normalement communiqués dans un courriel par le laboratoire du CQGC, mais peuvent aussi être récupérées dans la SampleSheet ou depuis BSSH (sous l'onglet "_Runs_"). Un modèle `README-${FC_SHORT}-template.ipynb` est disponible sur GitHub et dans Teams.

**_N.B._** Le format `ipynb` (Jupyter Notebook) est utilisé car du code peut y être exécuté, mais un simple format `txt` ou `md` peut aussi bien servir. Un fichier `README-FC_SHORT-template.ipynb` est disponible sur GitHub. Il contient le code à copier-coller et les cases à remplir pour le suivi des opérations.

#### 3. Se connecter à spxp-app02 

Au préalable, obtenir un jeton auprès du service informatique afin de pouvoir se connecter à distance _via_ le VPN du CHUSJ.

`ssh ${USER}@10.128.80.26`

#### 4. Mettre en place l'environnement de travail 

```bash
## 0. Mise en place de l'environnement de travail
# ssh ${USER}@10.128.80.26 # Renseigner la valeur de ${FC}

screen -S prag
conda activate CQGC-utils

export FC=""
a=($(echo ${FC} | tr '_' '\n'))
export FC_SHORT="${a[1]}_${a[2]}"
export BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/${a[1]}"
export WORKDIR="/mnt/spxp-app02/staging2/dragen"
```

### 1. Collecter les informations sur les familles et les métriques DragenGermline

Les informations sur la constitution des familles ("Case") sont centralisées dans Nanuq. La commande suivante permet de récupérer la liste des échantillons sur la _Run_ et de fournir les données nécessaires à la création de cas dans Emedgene, incluant les termes HPOs associées aux patients dans la base de données [Phenotips](https://chusj.phenotips.com//). Avec le NoveSeqX, il est également possible de générer les alignements localement et de récupérer les métriques pour faire approuver les échantillons avant que ceux-ci soient soumis aux analyses dans Emedgene.

***_N.B._***:
- Si le séquençage est réalisé avec le NovaSeq6000, il ne sera malheureusement pas possible de récupérer les métriques avant de créer les cas dans Emedgene. Dans ce cas, il faut s'assurer que le labo a demandé au NovaSeq6000 de faire la déconvolution sur BaseSpace et passer à l'étape 3, "créer les cas sur Emedgene", une fois que les FASTQs ont été générés (suivi sur l'onglet [Analyses](https://chusj.cac1.sh.basespace.illumina.com/analyses)). 
- Dans la situation où la déconvolution n'a pas lieu sur BaseSpace, il faut déconvoluer les BCLs avec les serveurs DRAGEN du CQGC et envoyer les FASTQs résultants sur BaseSpace (passer à l'étape 2, "Téléverser les FASTQs sur BaseSpace").

```bash
cd ${WORKDIR}
python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_collect_dragen_metrics.py ${FC}
```

_E.g._ `python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_collect_dragen_metrics.py 20240705_LH00336_0073_A22MFJFLT3`

La commande ci-dessus génère: 

- Les métriques des analyses DragenGermline en format HTML, `${FC_SHORT}_metrics.html`
- Les mêmes métriques en format CSV, `${FC_SHORT}_metrics.csv`
- Un fichier "samples_list.csv" qui sera utile pour les étapes suivantes. Voici un exmple du contenu de samples_list.csv

Le script s'attend à trouver les analyses dans le répertoire `/staging/hiseq_raw/LH00336/${FC}/Analysis/1/Data/DragenGermline/`.

Examiner les métriques et le fichier "samples_list.csv" générés par la commande ci-dessus. 

Envoyer les fichiers `${FC_SHORT}_metrics.csv` et `${FC_SHORT}_metrics.html` dans un courriel aux personnes responsables:

- "Dr Jacques Michaud (HSJ)" <jacques.michaud.med@ssss.gouv.qc.ca>
- "Camille Varin-Tremblay (HSJ)" <camille.varin-tremblay.hsj@ssss.gouv.qc.ca>
- "Rene Allard (HSJ)" <rene.allard.hsj@ssss.gouv.qc.ca>

Si les métriques passent les critères d'acceptabilité, procéder au téléversement des FASTQs dans BaseSpace. Sinon, attendre la réponse des responsables. Dans ce cas, il faudra éventuellement supprimer les lignes correspondantes aux échantillons de mauvaise qualité du fichier `samples_list.csv` avant de procéder aux prochaines étapes.

Seuils d'acceptabilité (en discussion):

- Coverage Uniformity < 0.25 ou 0.2
- Average Alignment Coverage > 25X ou 30X
- % mapped reads > 95%
- à partir des fichiers .wgs_coverage_metrics.csv:
    - COVERAGE SUMMARY,,Average alignment coverage over genome > 30X
    - COVERAGE SUMMARY,,PCT of genome with coverage [  20x: inf) > 90%
    - COVERAGE SUMMARY,,Uniformity of coverage (PCT > 0.4*mean) over genome > 90%
- À partir des fichiers .mapping_metrics.csv:
    - MAPPING/ALIGNING SUMMARY,,Number of unique & mapped reads (excl. duplicate marked reads) > 80%
    - MAPPING/ALIGNING SUMMARY,,Estimated sample contamination < 0.05
- À partir des fichiers .cnv_metrics.csv:
    - CNV SUMMARY,,Coverage uniformity
    - CNV SUMMARY,,Number of passing amplifications
    - CNV SUMMARY,,Number of passing deletions
- À partir des fichiers .ploidy_estimation_metrics.csv:
    - PLOIDY ESTIMATION,,Ploidy estimation doit correspondre au sexe de l'échantillon
- Autre chose qu'on devrait pouvoir confirmer avant de releaser les données est la structure familiale des trios

Exemple du contenu d'un fichier `samples_list.csv`:

    sample_name,biosample,relation,gender,ep_label,mrn,cohort_type,status,family_id,birthdate,flowcell_date,flowcell
    GM241567,27556,PROBAND,FEMALE,CHUSJ,03486257,TRIO,AFF,03486257,2024-04-29,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241601,27560,MTH,FEMALE,CHUSJ,03487612,TRIO,UNF,03486257,1980-10-15,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241575,27559,FTH,MALE,CHUSJ,03487451,TRIO,UNF,03486257,1978-10-02,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241566,27555,PROBAND,FEMALE,CHUSJ,03486541,TRIO,AFF,03486541,2024-02-13,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241572,27557,MTH,FEMALE,CHUSJ,03486957,TRIO,UNF,03486541,1987-05-04,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241573,27558,FTH,MALE,CHUSJ,03486959,TRIO,UNF,03486541,1987-01-20,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    24-05914-T1,27573,PROBAND,MALE,CHUS,2552542,TRIO,AFF,24-05914-T1,2024-06-09,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    24-05822-T1,27574,MTH,FEMALE,CHUS,24-05822-T1,TRIO,UNK,24-05914-T1,2000-11-28,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    24-05821-T1,27575,FTH,MALE,CHUS,24-05821-T1,TRIO,UNK,24-05914-T1,1989-01-24,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008162,27484,PROBAND,FEMALE,CUSM,MCH_5991033,TRIO,AFF,24-38625,2024-06-05,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008381,27485,MTH,FEMALE,CUSM,RVH_5994052,TRIO,UNF,24-38625,1986-10-13,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008380,27486,FTH,MALE,CUSM,R2265273,TRIO,UNF,24-38625,1994-06-01,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008229,27487,PROBAND,FEMALE,CUSM,MCH_5991739,TRIO,AFF,24-38627,2024-06-05,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008668,27488,MTH,FEMALE,CUSM,RVH_5175816,TRIO,UNF,24-38627,1994-09-29,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    MO-24-008667,27489,FTH,MALE,CUSM,RVH_5809722,TRIO,UNF,24-38627,1989-08-11,2024-07-05,20240705_LH00336_0073_A22MFJFLT3


### 2. Téléverser les FASTQs sur BaseSpace

À ce jour, il n'est pas possible de créer les cas dans Emedgene en partant de CRAM ou de VCF, tout en bénéficiant de toutes les fonctionnalités (ex: Manta, SMN callers,...). Il nous faut donc téléverser et repartir des FASTQs.

```bash
cd ${WORKDIR}/${FC_SHORT}
python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_upload_fastqs.py
```

Le script s'attend à trouver les fichiers dans le répertoire `/staging/hiseq_raw/LH00336/${FC}/Analysis/1/Data/DragenGermline/fastq`.


### 3. Créer les cas sur Emedgene

Cette étape consiste à:

1. Générer le fichier "emg_batch_manifest.csv"
2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
3. (**TODO**) Ajouter les participants _via_ l'API

#### 1. Générer le fichier "emg_batch_manifest.csv"

```bash
cd ${WORKDIR}
python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_make_batch_from_nanuq.py ${FC_SHORT} 2>&1 | tee ${FC_SHORT}/emg_make_batch.log
```

**_N.b._**: Les termes HPOs [Human Phenotype Ontology](https://hpo.jax.org/app/) sont contenus dans la base de données patients, Phenotips, et peuvent être consultées _via_ l'API avec un identifiant externe composé du "site + MRN" (ex. "CHUS1636084"). C'est ce que fait le script `emg_make_batch_from_nanuq.py` car l'identifiant primaire PhenotipsID (PID), n'est pas contenu dans Nanuq.

Alternativement, les informations sur les familles pour la _Run_, incluant le PID, peuvent être récupérées à partir du fichier Excel partagé sur "\\shsjcifs01\public crhsj\centre genomique\Projets\Pragmatiq\Étude PRAGMatIQ - Base de données.xlsx" ou sur Teams/OneDrive "Étude PRAGMatIQ - Base de données.xlsx" (canal de l'équipe _CHUSJ - Étude PRAGMatIQ_ dans Teams ("Documents > Général")). Ceci peut s'avérer nécessaire dans l'éventualité où le script `emg_make_batch_from_nanuq.py` échoue. _C.f._ Annexe A (Rechercher les termes HPOs) plus bas. En dernier recours, contacter Camille VARIN-TREMBLAY pour obtenir les informations manquantes.

#### 2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"

Le script Python `emg_make_batch_from_nanuq.py` exécuté à l'étape précédente génère le fichier manifeste (CSV) d'entrée "emg_batch_manifest.csv", avec les chemins d'accès complets aux FASTQs sur BSSH. Pour plus d'information au sujet du fichier d'entrée pour la création des cas par lot, consulter les spécifications d'Emedgene [CSV manifest specification](https://help.emedgene.com/en/articles/7231644-csv-format-requirements).

Il faut à présent se connecter _via_ l'interface web à Emedgene pour téléverser le fichier manifestedécrivant les cas à l'étude, puis lancer les analyses.

1. Se connecter à [Emedgene](https://chusaintejustine.emedgene.com/) avec le compte cqgc.bioinfo.hsj@ssss.gouv.qc.ca.
2. Cliquer sur "+ New Case", en haut de la page à droite
3. Cliquer sur "Switch to batch"
4. Glisser-déposer le fichier CSV manifest dans la boîte "Browse files from your computer or drag it here ('csv' | Maximum file size: 5 MB, 50 cases)"
5. Cliquer sur "Launch Case".

Pour plus d'informations, voir les instructions d'_Emedgene_ sur comment [créer des cas par lot](https://r4a56nl8uxkx3w3a292kjabk.emedgene.com/articles/7221986-batch-case-upload). https://chusaintejustine.emedgene.com/v2/#/help/

Exemple d'erreur fréquent: "Unmatched phenotype(s) found: phenotypes/EMG_PHENOTYPE_0011438". Il faut supprimer ce phénotype de la liste. Pour trouver le terme fautif, il faut remplacer "EMG_PHENOTYPE_" par "HP:" (dans ce cas-ci, EMG_PHENOTYPE_0011438 == HP:0011438)


#### 3. (**TODO**) Ajouter les participants


### 4. Archiver les résultats

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


### 5. Nettoyer

Une fois que les analyses sont terminées et archivées, il faut libérer l'espace disque sur `spxp-app02`. Au préalable, notifier le laboratoire, qui doit sauvegarder quelqeues métriques de l'opération de séquençage avant qU'elles ne soient effacées.


## Résolution des problèmes fréquents

### HPOs non-renseignés

Isoler les individus du cas concernés dans un fichier `emg_batch_manifest-[cas-sans-hpos].csv`. Une fois que le dossier Phenotips sera rempli, utiliser le script pour récupérer les termes HPOs à copier-coller dans le manifest. Puis, déposer le manifest dans Emedgene 

_C.f. Annexe A: Rechercher les termes HPOs_ 


## Références

- [Production environment](https://chusaintejustine.emedgene.com)
- [Eval environment](https://stejustine.emedgene.com/)
- [Phenotips](https://chusj.phenotips.com/)
- [Phenotips API](https://docs.phenotips.com/reference/accessandauthentication)
- [Human Phenotype Ontology (HPOs)](https://hpo.jax.org/app/)
- [Batch case upload instructions by Emedgene](https://r4a56nl8uxkx3w3a292kjabk.emedgene.com/articles/7221986-batch-case-upload) 
- https://chusaintejustine.emedgene.com/v2/#/help/
- [CSV manifest specification](https://help.emedgene.com/en/articles/7231644-csv-format-requirements)
- Old/Alternative method, using scripts:
    - "Case_creation_script_v2.docx"
    - `create_batch_cases_v2.py --help`
    - "create_batch_cases_v2_template.csv"


### Liste des scripts et leur usage

- `archive_PRAGMatIQ.sh`: Copier les résultats des analyses PRAGMatIQ du AWS d'Emedgene vers Narval.
- `create_batch_cases_v2.py`:  OBSOLETE - Script d'Emedgene pour soumettre des lots de cas (batch), utilisé à titre d'exemple de code.
- `emg_collect_dragen_metrics.py`: Générer les métriques pour les échantillons d'une _Run_ et créer le fichier `samples_list.csv`.
- `emg_collect_samples_metrics_all.py`: Générer les métriques pour _tous_ les échantillons dans le dossier d'archivage sur Narval.
- `emg_collect_samples_metrics.py`: Générer les métriques pour les échantillons d'une _Run_ à partir des résultats sur AWS.
- `emg_create_cases.py`: DEV - Prototype pour créer des cas à partir de JSON et des appels API à Emedgene.
- `emg_make_batch_from_nanuq.py`: Générer le manifeste pour créer des cas par lot dans Emedgene, à partir des donnes SampleNames de Nanuq.
- `emg_make_batch.py`: Générer le manifeste pour créer des cas par lot dans Emedgene, à partir d'une `samples_list.csv`.
- `emg_upload_fastqs.py`: Téléverser des FASTQs sur BaseSpace


## Annexe

Créer les cas avec la commande (ne fonctionne pas bien). Privilégier l'interface web "Switch to batch" pour y déposer le fichier manifeste. `python create_batch_cases_v2.py -i emg_batch_manifest.csv -s 10123 -hu stejustine.emedgene.com -u cqgc.bioinfo.hsj@ssss.gouv.qc.ca -p 7TmbuM3TUCMwP -b`.


### A. Rechercher les termes HPOs

Dans la situation où les PID doivent être récupérés manuellement pour chaque cas, utiliser le script `/staging2/soft/CQGC-utils/Helpers/get_phenotips_hpos.py`. Par exemple:

```bash
python /staging2/soft/CQGC-utils/Helpers/get_phenotips_hpos.py P0000844

# HP:0000403      Recurrent otitis media
# HP:0001935      Microcytic anemia
# HP:0004372      Reduced consciousness/confusion
# HP:0007204      Diffuse white matter abnormalities
# HP:0007305      CNS demyelination

# Format for Emedgene's WebUI:
# HP:0000403,HP:0001935,HP:0004372,HP:0007204,HP:0007305

# Format for Emedgene's batch:
# HP:0000403;HP:0001935;HP:0004372;HP:0007204;HP:0007305
```

Ou directement depuis un Notebook `ipynb`, copier-coller le code suivant dans une cellule

```python
import sys
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


### C. aws

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
