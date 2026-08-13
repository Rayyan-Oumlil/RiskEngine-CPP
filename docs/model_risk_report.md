# RiskEngine-CPP — Rapport de risque de méthode, de mesure et de modèle

> **Statut :** en cours de rédaction. §2 et §3 sont rédigés (phases 1 et 2). Les autres sections
> sont remplies au fil des phases du plan ([`riskengine_research.md`](riskengine_research.md) §12).
>
> **Règles éditoriales.** Chaque figure a une légende autonome (*ce qu'elle montre*, puis *pourquoi
> c'est important*). Chaque chiffre renvoie au CSV ou au test qui le produit. Chaque estimation
> stochastique est donnée avec son erreur standard. Chaque affirmation négative (« X échoue ») est
> accompagnée du mécanisme et de la parade.

## Résumé exécutif

*À rédiger en phase 8, une fois les découvertes chiffrées de §3 à §8 établies.*

---

## 1. Introduction et taxonomie du risque

*À rédiger en phase 8.* Cadre : [`riskengine_research.md`](riskengine_research.md) §0.

---

## 2. Cadre, conventions et reproductibilité

### 2.1 Dynamiques, paramètres de référence, conventions d'unités

Les unités (taux continus décimaux, vol annualisée décimale, maturités en années, Grecs bruts par
unité d'entrée) et le traitement des cas dégénérés sont fixés dans
[`conventions.md`](conventions.md). Le vecteur de référence de ce rapport est S = K = 100,
r = 5 %, q = 0, σ = 20 %, T = 1 an.

### 2.2 Protocole : générateur, graines, métrique d'efficacité, régénération

**Générateur.** Philox 4×32-10, un générateur à compteur : le tirage *i* du bloc *b* vaut
Philox(clé = graine, compteur = (i, b, flux)). Il n'y a aucun état à partager ni à faire avancer.
L'implémentation est vérifiée à la compilation contre les vecteurs de référence de Random123. Les
uniformes sont les milieux (k + ½)·2⁻⁵² d'une grille de 2⁵² points : jamais 0 ni 1, et 1 − u est
exact. Les normales sont obtenues par inversion (Wichura AS241, erreur relative ≤ 10⁻¹⁵ contre une
référence à 60 chiffres), jamais par `std::normal_distribution`, dont l'algorithme dépend de la
bibliothèque standard, ni par Box-Muller, qui détruit la structure des points quasi-aléatoires. La
grille étant symétrique, `normal(1 − u) == −normal(u)` exactement.

**Reproductibilité bit à bit.** Une simulation est découpée en un nombre fixe de blocs, indépendant
du nombre de threads. Chaque bloc possède son flux aléatoire et son accumulateur de Welford ; les
résultats partiels sont fusionnés (Chan) dans l'ordre des blocs. L'addition flottante n'étant pas
associative, c'est cet ordre fixe qui garantit un résultat identique au bit près avec 1, 2, 3, 8 ou
64 threads. La contraction automatique en FMA est désactivée (`-ffp-contract=off`), si bien que GCC
et Clang produisent les mêmes bits, y compris avec `-march=native`. Une estimation Monte Carlo est
donc une fonction pure de (graine, nombre de tirages, nombre de blocs).

**Portes de validation (phase 2)**, toutes dans `tests/test_rng.cpp` et `tests/test_simulation.cpp`
avec graines figées (tests déterministes) :

| Porte | Résultat |
|---|---|
| Philox contre les vecteurs Random123 | 3/3, vérifié par `static_assert` |
| Uniformes : moyenne, variance, Kolmogorov-Smirnov (10⁶ tirages) | dans 4 écarts-types ; √n·D < 1,95 (seuil 99,9 %) |
| Normales : moyenne, variance, asymétrie, kurtosis, Kolmogorov-Smirnov (10⁶ tirages) | idem |
| Corrélation sérielle (retard 1) et entre blocs | < 4/√n |
| Résultat identique à 1/2/3/8/64 threads | égalité exacte |
| Même résultat sur GCC et Clang | égalité exacte sur valeurs de référence (golden) |
| Calibration de l'erreur standard, 200 réplications indépendantes | ≥ 90 % des \|z\| ≤ 1,96, aucun \|z\| > 4 |

**Efficacité.** Les techniques de réduction de variance seront comparées par
efficacité = 1 / (variance × temps CPU) (Glasserman), à partir de la phase 4.

**Régénération.** Chaque figure ou table est produite par un exécutable de `experiments/`, qui écrit
`data/results/<id>.csv` et `<id>.meta.json` : commit git et indicateur de copie de travail modifiée,
compilateur, configuration, drapeaux, paramètres, date. Les résultats versionnés proviennent d'un
build Release d'un arbre propre (`"git_dirty": false`). Les figures sont produites par
`tools/make_figures.py` :

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
build-release/experiments/fd_vcurve && build-release/experiments/iv_roundtrip
python3 tools/make_figures.py        # pip install -r tools/requirements.txt
```

---

## 3. Vérité terrain : Black-Scholes analytique

Toutes les méthodes des sections suivantes (arbres, Monte Carlo, Grecs stochastiques) sont jugées
contre la formule fermée. Cette section établit que la formule fermée elle-même est exacte à la
précision affichée, et documente les deux limites numériques qu'elle impose déjà :
l'inversion en volatilité implicite, et les différences finies.

### 3.1 Validation, invariants, volatilité implicite

**Valeurs de référence.** Dix cas (vecteur canonique, dividende continu, exemple 15.6 de Hull, et
deux strikes d'aile dont le prix hors de la monnaie vaut ~10⁻¹²) sont comparés à une référence mpmath
à 50 chiffres ([`tools/bs_reference.py`](../tools/bs_reference.py)). Les Grecs de référence y sont
des *dérivées numériques du prix* à haute précision, et non les formules fermées : ils vérifient donc
indépendamment les formules, leurs signes et le traitement du dividende. Le prix et les cinq Grecs
concordent à 10⁻⁹ en relatif dans les dix cas, y compris pour les prix de 10⁻¹² en aile
(`tests/test_black_scholes.cpp`).

| Vecteur canonique | Valeur |
|---|---|
| Call / put | 10,450584 / 5,573526 |
| Delta call | 0,636831 |
| Gamma | 0,018762 |
| Vega | 37,524 par 1,00 de vol (0,375 par point) |
| Theta call | −6,414 par an (−0,01757 par jour calendaire) |
| Rho call | 53,232 par 1,00 de taux |

**Invariants**, sur une grille de 315 points (strike 50 à 200 pour S = 100, maturité 1 jour à 5 ans,
vol 5 % à 100 %) : parité call-put à 10⁻¹⁴ près relativement à S + K ; bornes de non-arbitrage ;
∂C/∂K < 0 et convexité en strike ; Grecs fermés conformes aux différences finies du prix. Pour
T → 0 et σ → 0, le pricer renvoie la valeur intrinsèque actualisée du forward par une branche
explicite, jamais par une division par zéro ; aucun NaN n'apparaît jusqu'à T = σ = 10⁻³⁰⁰.

**Volatilité implicite.** Le solveur ramène la cotation à l'option hors de la monnaie par la parité,
puis résout log prix(σ) = log cible par la méthode de Brent sur un intervalle encadrant. Travailler
en log prix garde le problème bien conditionné dans les ailes, où le prix couvre des centaines
d'ordres de grandeur. Newton seul n'est pas utilisé : il diverge là où vega → 0.

![Erreur relative de la vol implicite retrouvée en fonction de la borne de bruit de la cotation](figures/iv_roundtrip.svg)

*Figure 1 — Aller-retour prix → vol → prix sur la grille de 315 points, pour la cotation hors de
la monnaie (bleu) et dans la monnaie (orange) de chaque point. En abscisse, la borne de bruit
ε·(1 + (1 + d²)(a + b)/(σ·vega)), où a − b est la formule du prix : c'est l'erreur relative de
vol qu'entraîne à elle seule l'erreur d'arrondi de la cotation. **Ce qu'elle montre :** tous les
points sont sous la diagonale ; l'erreur du solveur n'excède jamais le bruit de la cotation, et la
dégradation dans la monnaie suit exactement ce bruit. **Pourquoi c'est important :** une vol
implicite imprécise tirée d'une cotation dans la monnaie n'est pas un défaut du solveur mais une
perte d'information dans la cotation elle-même ; aucun solveur ne peut la récupérer.*
Source : [`data/results/iv_roundtrip.csv`](../data/results/iv_roundtrip.csv).

| Cotation | Résolus | Erreur relative max | Au-delà de 10⁻⁹ | Erreur / borne, max |
|---|---|---|---|---|
| Hors de la monnaie | 289 / 315 | 1,3 × 10⁻¹³ | 0 | 0,68 |
| Dans la monnaie | 226 / 315 | 5,7 × 10⁻² | 19 | 0,64 |

Les 26 cotations hors de la monnaie non résolues ont un prix nul par dépassement de capacité
inférieur (par exemple 1 jour, K/S = 2, σ = 5 %) : il n'y a pas de vol à retrouver. Les 89
cotations dans la monnaie non résolues ont une valeur temps inférieure à l'arrondi de leur valeur
intrinsèque ; le solveur le signale (`ZeroTimeValue`) au lieu de renvoyer un nombre. Le pire cas
résolu, un put K = 150 à 3 mois et σ = 10 %, retrouve la vol à 5,7 % près.

**Mécanisme.** Un prix Black-Scholes est une différence a − b de deux termes (pour un call,
a = S e^{−qT} N(d₁) et b = K e^{−rT} N(d₂)). Dans la queue gaussienne, l'arrondi de d est amplifié
d'un facteur d, si bien que chaque terme porte une erreur relative d'environ ε·(1 + d²). La
cotation n'est donc connue qu'à ε·(1 + d²)·(a + b) près. Hors de la monnaie, a + b reste du même
ordre que le prix et la vol est retrouvée à ~10⁻¹³. Dans la monnaie, a + b ≈ S + K alors que la
valeur temps est minuscule : la perte atteint plusieurs pour cent.

**Parade.** Toujours inverser sur la cotation hors de la monnaie ; une vol implicite tirée d'une
cotation profondément dans la monnaie doit être accompagnée de sa borne de bruit. La formulation de
Jäckel (*Let's Be Rational*, 2015) reste la cible pour réduire la part due au pricer lui-même dans
les ailes.

### 3.2 Arrondi flottant et courbe en V des différences finies

![Erreur absolue du delta et du gamma par différences finies en fonction du bump relatif](figures/fd_vcurve.svg)

*Figure 2 — Erreur absolue du delta et du gamma par différences centrées sur le prix Black-Scholes
(vecteur canonique), contre la formule fermée, pour un bump relatif h de 10⁻¹⁴ à 10⁻¹. Les
pointillés marquent les optimums théoriques ε^{1/3} (delta) et ε^{1/4} (gamma). **Ce qu'elle
montre :** l'erreur a la forme d'un V ; à droite, la troncature décroît en h² ; à gauche, l'arrondi du
prix croît en ε/h pour le delta et en ε/h² pour le gamma. **Pourquoi c'est important :** diminuer le
bump ne rend pas un Grec plus précis ; en deçà de l'optimum, il le rend faux, et pour le gamma
catastrophiquement.* Source : [`data/results/fd_vcurve.csv`](../data/results/fd_vcurve.csv).

| Bump relatif h | Erreur delta | Erreur gamma (relative) |
|---|---|---|
| 10⁻³ | 8,6 × 10⁻⁷ | 1,2 × 10⁻⁶ |
| 10⁻⁴ | 8,6 × 10⁻⁹ | 1,4 × 10⁻⁸ |
| 10⁻⁵ | 8,1 × 10⁻¹¹ | 7,5 × 10⁻⁸ |
| 10⁻⁸ | 1,7 × 10⁻⁹ | 24 % |
| 10⁻¹⁴ | 1,1 × 10⁻³ | 7,7 × 10¹¹ |

Le plancher d'erreur du delta est d'environ 10⁻¹¹ pour h entre 3 × 10⁻⁷ et ε^{1/3} ≈ 6 × 10⁻⁶ ;
celui du gamma d'environ 10⁻¹⁰ (10⁻⁸ en relatif) autour de h ∈ [10⁻⁵, 10⁻⁴], près de
ε^{1/4} ≈ 1,2 × 10⁻⁴. Les pics isolés vers le bas (par exemple 4 × 10⁻¹⁵ pour le delta à
h = 2,4 × 10⁻⁶) sont des changements de signe de l'erreur, pas une précision atteignable. En
deçà de h ≈ 5 × 10⁻⁹, les 46 points de la grille donnent un gamma faux à 100 % ou plus, dont 15
exactement nuls : le numérateur p(S + h) − 2p(S) + p(S − h) s'annule au bit près.

**Mécanisme.** Pour la différence centrée, l'erreur de troncature vaut ~(h²S²/6)·∂³V/∂S³ et
l'erreur d'arrondi ~ε·V/(hS) pour le delta, ~ε·V/(hS)² pour le gamma ; leur somme est minimale
en h* ∝ ε^{1/3} et ε^{1/4} respectivement. L'expérience divise par le pas effectivement réalisé en
flottant, et non par le pas voulu, afin que la courbe ne montre que l'arrondi du prix.

**Parade.** Bump relatif, jamais absolu, et proche de l'optimum de l'ordre de dérivation. Un bump
de 10⁻⁴ en relatif, choix courant, est sûr pour les deux Grecs sur un pricer déterministe. Sur un
pricer Monte Carlo, le bruit statistique remplace l'arrondi et déplace l'optimum de plusieurs ordres
de grandeur : c'est l'objet de §6.

---

## 4. Arbres : convergence et pathologies

*À rédiger en phase 3.*

## 5. Monte Carlo : convergence et réduction de variance

*À rédiger en phase 4.*

## 6. Estimation des Grecs sous bruit *(section phare)*

*À rédiger en phase 5.*

## 7. Mesures de risque sur positions non linéaires

*À rédiger en phase 6.*

## 8. Risque de modèle au-delà de GBM

*Phase 7, optionnelle ([`riskengine_research.md`](riskengine_research.md) §12.1).*

## 9. Notes d'ingénierie et performance

*À rédiger en phase 8.*

## 10. Limites et travaux futurs

*À rédiger en phase 8.*

## 11. Conclusion : recommandations opérationnelles

*À rédiger en phase 8.*

---

## Annexes

*A. dérivations des estimateurs pathwise et LR ; B. tables de paramètres ; C. valeurs de référence
et sources ; D. procédure de reproduction (voir §2.2 pour l'état actuel).*
