# Procédure PRAGMatIQ

Cette procédure décrit comment créer des cas (trio, solo, duo ou quad) à analyser sur la plateforme d'[Emedgene (EMG)](https://chusaintejustine.emedgene.com/) une fois que les séquences des échantillons sont disponibles sur BaseSpace Sequence Hub (BSSH). Les BCLs issus des séquenceurs sont d'abord convertis en FASTQs par DRAGEN BCLConvert, soit directement sur le NovaSeqX, soit par notre serveur local ou sur BSSH. La progression des opérations de séquençage puis de la déconvolution-conversion de la _Run_ peuvent être suivies sur BSSH sous l'onglet ["Run"](https://chusj.cac1.sh.basespace.illumina.com/runs/active) et ["Analyses"](https://chusj.cac1.sh.basespace.illumina.com/analyses/), respectivement. Les cas à analyser sont finalement créés dans Emedgene à l'aide d'un fichier "batch manifest" en format CSV. Ce fichier contient les informations sur les individus des familles (sexe, âge, relation, FASTQs associés, termes HPOs, _etc_.). 

En résumé, voici les étapes de la création des cas sur Emedgene pour une run donnée ("flowcell", `${FC}`). La déconvolution-conversion des BCLs en FASTQs (étape 1) est normalement prise en charge par un processus en `cron`, `dragen_bcl-convert_watcher.sh`. Le script `run_ipeline_prag.sh` réalise quant à lui les étapes 2 à 4.1.

0. Pré-requis
    1. Obtenir les identifiants de connexion
    2. Créer un journal pour le suivi, par exemple `README-${FC_SHORT}.ipynb` (*)
    3. Se connecter à spxp-app02 `ssh ${USER}@10.128.80.26`
    4. Mettre en place l'environnement de travail `conda activate CQGC-utils`.
1. Déconvolution et conversion des BCLs en FASTQs `dragen_bcl-convert_launcher.sh`
2. Collecter les informations sur les familles `list_run_samples.py ${FC}`.
3. Téléverser les FASTQs sur BaseSpace `emg_upload_fastqs.py`
4. Créer les cas sur Emedgene 
    1. Générer le fichier "emg_batch_manifest.csv" `emg_make_batch.py`
    2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
    3. (**TODO**) Ajouter les participants _via_ l'API
5. Collecter les métriques `emg_collect_samples_metrics.py`
6. Archiver les résultats
7. Nettoyer

Où `${FC_SHORT}` est le nom court de la _FlowCell/Run_. Ex: Si la _flowcell/Run_ se nomme "230727_A00516_0441_AHKVFYDMXY", `${FC_SHORT}` est "A00516_0441".

(*) Un fichier `README-FC_SHORT-template.ipynb` est disponible sur [GitHub](https://github.com/CQGC-Ste-Justine/CQGC-utils/tree/main/Analysis.pipeline_illumina). Il contient le code à copier-coller et les cases à remplir pour le suivi des opérations.

***_N.B._***: Afin de connecter les informations génétiques des familles (issues de Nanuq) aux phénotypes (termes HPOs contenus dans Phenotips), il est impératif que les deux champs EP (Établissement Public) et MRN (Medical Record Number) soient bien renseignés dans les deux systèmes par les collègues. Sinon, il faut obtenir les informations de Phenotips (identifiants PID) par "Camille Varin-Tremblay (HSJ)" <camille.varin-tremblay.hsj@ssss.gouv.qc.ca>.


## Procédure


### 0. Pré-requis

Installation des outils et obtention des identifiants de connexion aux différents services.

#### 0.1. Obtenir les identifiants de connexion

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
        (...)
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

#### 0.2. Créer un journal pour le suivi

Créer un journal dans lequel seront notés le suivi des opérations. Par exemple, dans un notebook Jupyter nommé `README-${FC_SHORT}.ipynb` et y inscrire en titre les noms de la _flowcell_ (_e.g._.: "230711_A00516_0433_BHKVMJDMXY", où `${FC_SHORT}` serait dans ce cas "A00516_0433") et de l'exérience ("Seq_S2_PRAG_20230711"). Ces informations sont normalement communiqués dans un courriel par le laboratoire du CQGC, mais peuvent aussi être récupérées dans la SampleSheet ou depuis BSSH (sous l'onglet "_Runs_"). Un modèle `README-${FC_SHORT}-template.ipynb` est disponible sur [GitHub](https://github.com/CQGC-Ste-Justine/PrivateDoc/blob/0f59d674fd5d91ca5da7880ab55cd753ad324203/gapp_conf.json) et dans Teams.

```bash
# Lister le nom de la Flowcell (FC) et de l'expérience (XP)
bs -c cac1 list runs --newer-than 3d --filter-term PRAG --filter-field ExperimentName
```

**_N.B._** Le format `ipynb` (Jupyter Notebook) est utilisé car du code peut y être exécuté, mais un simple format `txt` ou `md` peut aussi bien servir. Un fichier `README-FC_SHORT-template.ipynb` est disponible sur [GitHub](https://github.com/CQGC-Ste-Justine/PrivateDoc/blob/0f59d674fd5d91ca5da7880ab55cd753ad324203/gapp_conf.json). Il contient le code à copier-coller et les cases à remplir pour le suivi des opérations.

#### 0.3. Se connecter à spxp-app02 

Au préalable, obtenir un jeton auprès du service informatique afin de pouvoir se connecter à distance _via_ le VPN du CHUSJ.

`ssh ${USER}@10.128.80.26`

#### 0.4. Mettre en place l'environnement de travail 

```bash
## 0. Mise en place de l'environnement de travail
# ssh ${USER}@10.128.80.26
# Pour connaître l'identifiant FC et XP
bs -c cac1 list runs --newer-than 3d --filter-term PRAG --filter-field ExperimentName
screen -S prag

conda deactivate && conda activate CQGC-utils
export FC="" # Ex.: 20250613_LH00336_0223_A232CFLLT3
a=($(echo ${FC} | tr '_' '\n')) 
export BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/${a[1]}"
export WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs"
export SOFTDIR='/mnt/spxp-app02/staging2/soft/CQGC-utils'
```

### 1. Déconvolution et conversion des BCLs en FASTQs (automatique)

***N.B.*** Cette étape est normalement réalisée de manière automatique par le script `dragen_bcl-convert_watcher.sh` en `cron`. Le watcher surveille et vérifie certains fichiers marqueurs dans le dossier des BCLs (sorties des séquenceurs) avant de lancer `dragen_bcl-convert_launcher.sh`

Le script `run_pipeline_prag.sh`, qui est exécuté ici, surveille le répertoire des BCLs pour les fichiers marqueurs `CopyComplete.txt` (fin du séquençage) et `FastqCoplete.txt` (complétion de la déconvolution-conversion). 

```bash
if [[ ! -d "${WORKDIR}/${FC}" ]]; then mkdir ${WORKDIR}/${FC}; fi
cd ${WORKDIR}/${FC}
bash ${SOFTDIR}/Analysis.pipeline_illumina/run_pipeline_prag.sh ${FC} 2>&1 | tee ${WORKDIR}/${FC}/run_pipeline_prag.log
```

Les étapes 1 à 4.1 de cette procédure sont incluses dans le pipeline. Elles sont décrites ci-bas à titre informatif. À la sortie du pipeline, un fichier `emg_batch_manifets.csv` est généré. Ce document, qui est à glisser-déposer sur le site d'Emedgene, permet de créer les cas à analyser (passez à l'étape 4.2)

### 2. Collecter les informations sur les familles (`run_pieplie_prag.sh`)

Générer le fichier `samples_list.csv`.

```bash
## 2. Collecter les informations sur les familles dans samples_list.csv,
## fournit aussi la liste des échantillons à téléverser sur BaseSpace (étape 3)
echo "Get list of samples for run ${FC}"
if [[ ! -d "${WORKDIR}/${FC}" ]]; then
    mkdir ${WORKDIR}/${FC};
fi
cd ${WORKDIR}/${FC}
python ${SOFTDIR}/Analysis.pipeline_illumina/list_run_samples.py ${FC}
```

Les informations sur la constitution des familles ("Case") sont centralisées dans Nanuq. La commande ci-dessus permet de récupérer la liste des échantillons sur la _Run_ et de fournir les données nécessaires à la création de cas dans Emedgene, incluant les termes HPOs associées aux patients dans la base de données [Phenotips](https://chusj.phenotips.com//). 

Avec le NoveSeqX, il est également possible de générer les alignements localement et de récupérer les métriques pour faire approuver les échantillons avant que ceux-ci soient soumis aux analyses dans Emedgene, mais cette solutio a été rejetée car les analyses exécutées directement sur le NovaSeqX bloqent le séquenceur qui ne peut pas séquencer d'autres flowcells.

Le fichier "samples_list.csv" sera utile pour les étapes suivantes. Voici un exmple du contenu de `samples_list.csv`:

    sample_name,biosample,relation,gender,ep_label,mrn,status,family_id,birthdate,project,flowcell
    25-05853-T1,36712,PROBAND,FEMALE,CHUQ,1778336,AFF,25-05853-T1,2023-07-21,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05855-T1,36713,MTH,FEMALE,CHUQ,782792,UNK,25-05853-T1,1988-03-29,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05858-T1,36714,FTH,MALE,CHUQ,1276400,UNK,25-05853-T1,1987-07-15,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05895-T1,36717,PROBAND,MALE,CHUS,2487041,AFF,25-05895-T1,2022-05-08,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05875-T1,36716,MTH,FEMALE,CHUS,2193537,UNK,25-05895-T1,2003-12-05,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05873-T1,36715,FTH,MALE,CHUS,2083226,UNK,25-05895-T1,2003-02-13,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-06003-T1,36711,PROBAND,FEMALE,CHUS,2591577,AFF,25-06003-T1,2025-05-24,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05903-T1,36709,MTH,FEMALE,CHUS,1497199,UNK,25-06003-T1,1992-07-10,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-05905-T1,36710,FTH,MALE,CHUS,25-05905-T1,UNK,25-06003-T1,1992-10-05,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-06095-T1,36884,PROBAND,MALE,CHUS,2591918,AFF,25-06095-T1,2025-05-28,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-06078-T1,36883,MTH,FEMALE,CHUS,635210,UNK,25-06095-T1,1993-11-15,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3
    25-06077-T1,36882,FTH,MALE,CHUS,371989,UNK,25-06095-T1,1992-04-06,PRAGMATIQ_CHUS,20250613_LH00336_0223_A232CFLLT3


### 3. Téléverser les FASTQs sur BaseSpace

ATTENTION! Si cette étape échoue ou est interrompue, il faut effacer les fichiers créés sur BaseSpace avant de relancer la commande. Pour effacer les fichiers FASTQs des biosamples sur BaseSpace, utilisez [l'interface web](https://chusj.cac1.sh.basespace.illumina.com/biosamples).

```bash
## 3. Téléverser les FASTQs sur BaseSpace
echo "Uploading samples to BaseSpace"
## TODO: move FASTQ files from 1.fastq/PROJECT_NAME to 1.fastq/ before uploading when instrument is NovaSeq6000
python ${SOFTDIR}/Analysis.pipeline_illumina/emg_upload_fastqs.py
touch ${WORKDIR}/${FC}/UploadBsComplete.txt
```

Il n'est pas possible de créer les cas dans Emedgene en partant de CRAM ou de VCF, tout en bénéficiant de toutes les fonctionnalités (ex: Manta, SMN callers,...). Il nous faut donc téléverser et repartir des FASTQs.

`emg_upload_fastqs.py` prend comme entrée le fichier `samples_list.csv`, créé à l'étape précédente.Sinon, il faut spécifier le chemin d'accès à un samples_list avec l'option `--file`.

Le script s'attend à trouver les fichiers dans le répertoire `${WORKDIR}/${FC}/1.fastq`, sauf si spécifié autrement avec l'option `--data-dir`.


### 4. Créer les cas sur Emedgene

Cette étape consiste à:

1. Générer le fichier "emg_batch_manifest.csv"
2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
3. (**TODO**) Ajouter les participants _via_ l'API

#### 4.1. Générer le fichier "emg_batch_manifest.csv"

```bash
## 4. Créer les cas sur Emedgene 
###  4.1. Générer le fichier "emg_batch_manifest.csv" `emg_make_batch_from_nanuq.py ${FC}`
###  4.2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
# TODO: use `bs` to check if files are on BaseSpace for each sample
# Fails if launched immediately after UploadBsComplete.
# Wait a while for the last upload to register with BaseSpace.
sleep ${NAPTIME}
python ${SOFTDIR}/Analysis.pipeline_illumina/emg_make_batch.py >> ${WORKDIR}/${FC}/emg_make_batch.log 2>&1
```

**_N.b._**: Les termes HPOs [Human Phenotype Ontology](https://hpo.jax.org/app/) sont contenus dans la base de données patients, Phenotips, et peuvent être consultées _via_ l'API avec un identifiant externe composé de l'"Établissement Public + MRN" (ex. "CHUS1636084"). C'est ce que fait le script `emg_make_batch_from_nanuq.py` car l'identifiant primaire PhenotipsID (PID), n'est pas contenu dans Nanuq.

Alternativement, les informations sur les familles pour la _Run_, incluant le PID, peuvent être récupérées à partir du fichier Excel partagé sur "\\shsjcifs01\public crhsj\centre genomique\Projets\Pragmatiq\Étude PRAGMatIQ - Base de données.xlsx" ou sur Teams/OneDrive "Étude PRAGMatIQ - Base de données.xlsx" (canal de l'équipe _CHUSJ - Étude PRAGMatIQ_ dans Teams ("Documents > Général")). Ceci peut s'avérer nécessaire dans l'éventualité où le script `emg_make_batch_from_nanuq.py` échoue. _C.f._ Annexe A (Rechercher les termes HPOs) plus bas. En dernier recours, contacter Camille VARIN-TREMBLAY pour obtenir les informations manquantes.

#### 4.2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"

Le script Python `emg_make_batch_from_nanuq.py` exécuté à l'étape précédente génère le fichier manifeste (CSV) d'entrée "emg_batch_manifest.csv", avec les chemins d'accès complets aux FASTQs sur BSSH. Pour plus d'information au sujet du fichier d'entrée pour la création des cas par lot, consulter les spécifications d'Emedgene [CSV manifest specification](https://help.emedgene.com/en/articles/7231644-csv-format-requirements).

Il faut à présent se connecter _via_ l'interface web à Emedgene pour téléverser le fichier manifestedécrivant les cas à l'étude, puis lancer les analyses.

1. Se connecter à [Emedgene](https://chusaintejustine.emedgene.com/) avec le compte cqgc.bioinfo.hsj@ssss.gouv.qc.ca.
2. Cliquer sur "+ New Case", en haut de la page à droite
3. Cliquer sur "Switch to batch"
4. Glisser-déposer le fichier CSV manifest dans la boîte "Browse files from your computer or drag it here ('csv' | Maximum file size: 5 MB, 50 cases)"
5. Cliquer sur "Launch Case".

Pour plus d'informations, voir les instructions d'_Emedgene_ sur comment [créer des cas par lot](https://r4a56nl8uxkx3w3a292kjabk.emedgene.com/articles/7221986-batch-case-upload). https://chusaintejustine.emedgene.com/v2/#/help/

Si Emedgene détecte des erreurs, télécharger le fichier CSV de validation, puis corriger les erreurs qui y sont énummérés dans la première colonne du CSV. Voir la section "Résolution des problèmes" plus bas.

#### 4.3. Ajouter les participants

(**TODO**) Ajouter les participants


### 5. Collecter les metriques

```bash
5. Collecter les metriques
python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_collect_samples_metrics.py ${FC}
```

Envoyer les fichiers `${FC}_metrics.csv` et `${FC}_metrics.html` dans un courriel aux personnes responsables:

- "Dr Jacques Michaud (HSJ)" <jacques.michaud.med@ssss.gouv.qc.ca>
- "Camille Varin-Tremblay (HSJ)" <camille.varin-tremblay.hsj@ssss.gouv.qc.ca>
- "Rene Allard (HSJ)" <rene.allard.hsj@ssss.gouv.qc.ca>

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

### 6. Archiver les résultats

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

### 7. Nettoyer

Une fois que les analyses sont terminées et archivées, il faut libérer l'espace disque sur `spxp-app02`. Au préalable, notifier le laboratoire, qui doit sauvegarder quelqeues métriques de l'opération de séquençage avant qU'elles ne soient effacées.


## Résolution des problèmes fréquents

Exemple d'erreur fréquents.

### "Unmatched phenotype(s) found"

"Unmatched phenotype(s) found: phenotypes/EMG_PHENOTYPE_0011438" 

Il faut supprimer ce phénotype de la liste. Pour trouver le terme fautif, il faut remplacer "EMG_PHENOTYPE_" par "HP:" (dans ce cas-ci, EMG_PHENOTYPE_0011438 == HP:0011438). Les termes HPOs associés à la mère ne sont pas acceptés par Emedgene et génèrent cet erreur (ex: [Maternal diabetes (HP:0009800)](http://purl.obolibrary.org/obo/HP_0009800)). 

### HPOs non-renseignés

Nous nous basons sur les informations du site (EP, établissement public) + le MRN (medical record number) pour faire le lien entre la base Nanuq (échantillons) et Phenotips (patients). Si EP ou MRN est mal renseigné dans l'une ou l'autre de ces bases d'informations, les termes HPOs ne peuvent pas être associés au cas. Le Helper script `get_phenotips_hpos.py` peut être utilisé pour recueillir les termes HPOs avec un PhenotipsID (PID) à posteriori. _C.f. Annexe A: Rechercher les termes HPOs_ 

Les termes HPOs peuvent être absents si:

- Erreur de saisie EP ou MRN dans Nanuq ou Phenotips. 
- Could not find HPO terms for PID=P0000844 (EP+MRN=CHUS2487041): Le dossier existe, mais le médecin n'a pas rempli les phénotypes.

Il faut contacter Camille Varin-Tremblay. En attendant, isoler les individus du cas concernés dans un fichier `emg_batch_manifest-[cas-sans-hpos].csv`. Une fois que le dossier Phenotips sera rempli, utiliser le script pour récupérer les termes HPOs à copier-coller dans le manifest. Puis, déposer le manifest dans Emedgene 

### FamilyID

Un membre sain de la famille (ex: père, mère) n'est pas connecté à un proband. Demander à Camille ou consulter la base de données Excel.

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
