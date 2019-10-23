# hyliter
R package to connect HyLiTE output to the DeSeq2 differential gene expression software

## Warning

this is extremely in-devlopment code, that may well not do what you want it to,
break unexpectedly or change in the future. If you still want to play with the
existing functions, you can instal the package using devtools

```{r}
devtools::install_github("dwinter/hyliter")
```

And use it, here H2 and H1 are hybirds with hylite output ditrectories and 

```{r}
parental_mod <- parental_DE ( read_hylite ("../H1.hylite/H1.hylite.expression.txt"), "P1", "P2")
H1_mod <- hybrid_DE ("../H1.hylite/", "H1")
H2_mod <- hybrid_DE ("../H2.hylite/", "H2")
```
