FROM python:3.11-slim

# Dependances systeme necessaires a la compilation de certains packages scientifiques
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Installation de inmoose (pyComBat-seq) et pandas
RUN pip install --no-cache-dir inmoose pandas

# Embarquer le script de correction dans l'image, HORS de /data :
# /data est le point de montage du volume utilisateur a l'execution,
# son contenu masquerait tout fichier COPY dans /data.
WORKDIR /app
COPY batch_correction.py .

# Repertoire de travail par defaut : le dossier des comptages monte en volume
WORKDIR /data

# Commande par defaut : le script de correction. Pour passer des arguments,
# il faut reecrire la commande complete apres le nom de l'image.
CMD ["python", "/app/batch_correction.py"]