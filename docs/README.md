Required software installed as follows:

```shell
$ mamba create -n docs python=3.12 sphinx==8.2.3 pydata-sphinx-theme==0.16.1 sphinx-design==0.6.1 sphinx-copybutton==0.5.2
$ mamba activate docs
```

To test document build:

```shell
$ sphinx-build -M html docs/ .test-docs/
```