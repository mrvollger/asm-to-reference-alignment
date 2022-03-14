# Reference alignment workflow

## an example run script

```
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --cores $threads --use-conda -p
```

Or if you want to distribute over a cluster:

```
mkdir dir -p logs/drmaa
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --jobs $threads --use-conda -p  \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree=8G -pe serial {threads} -V -cwd -S /bin/bash -w n" --drmaa-log-dir logs/drmaa
```

Or if you want to make ideograms:

```
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --cores $threads --use-conda -p ideogram
```

```mermaid
flowchart LR;
    A[THIS IS A];
    B(THIS IS B);
    C(THIS IS C);
    D(THIS IS D);
    d(THIS IS d);

    subgraph COMBO [name for combinations]
        direction LR;
        A-->D;
        A-->C;
        B-->D;
        B-->C;
    end

    start --> COMBO;
    COMBO --> END;
```

```mermaid
journey
    title My working day
    section Go to work
      Make tea: Me
      Go upstairs: 3: Me
      Do work: 1: Me, Cat
    section Go home
      Go downstairs: 4: Me
      Sit down: 5: Me
```

```mermaid
gantt
    title A Gantt Diagram
    dateFormat  YYYY-MM-DD
    section Section
        A task           :a1, 2014-1-1, 500d
        A task2           :after a1, 300d
```
