# RiskEngine-CPP — Plan de recherche et de développement (v2)

> **Statut :** plan unique et faisant autorité pour ce projet. Remplace v1.
> **Livrable unique :** `docs/model_risk_report.md`. Le code C++20 n'existe que pour produire, de façon reproductible, les preuves numériques de ce rapport.
> **Cadrage :** projet étudiant, construit à temps partiel à côté des cours. Le chemin critique (§12) est prioritaire sur toute extension ; la phase 7 (Heston/Merton) est explicitement la première chose à couper si le temps manque — voir §12.1.

**Note de provenance.** Les valeurs de référence marquées ✅ ont été recalculées indépendamment (Python/SciPy) lors de la rédaction de cette v2 et concordent avec la v1. Les chiffres de performance de la v1 (débits, speedups, false sharing) sont conservés comme *hypothèses à re-mesurer* : ils dépendent de la machine et seront régénérés par `bench/` avec la configuration CPU consignée. Les sources de données de la section 9 sont à revérifier (quotas, conditions d'usage) avant la phase 6.

---

## Table des matières

0. [Thèse et taxonomie du risque](#0-thèse-et-taxonomie-du-risque)
1. [Ce qui change par rapport à la v1](#1-ce-qui-change-par-rapport-à-la-v1)
2. [Architecture C++20](#2-architecture-c20)
3. [Fondations mathématiques et valeurs de référence](#3-fondations-mathématiques-et-valeurs-de-référence)
4. [Monte Carlo : échantillonnage et réduction de variance](#4-monte-carlo--échantillonnage-et-réduction-de-variance)
5. [Estimation des Grecs sous bruit (section phare)](#5-estimation-des-grecs-sous-bruit-section-phare)
6. [Mesures de risque sur positions non linéaires](#6-mesures-de-risque-sur-positions-non-linéaires)
7. [Risque de modèle au-delà de GBM](#7-risque-de-modèle-au-delà-de-gbm)
8. [Ingénierie de performance](#8-ingénierie-de-performance)
9. [Données historiques et stress tests](#9-données-historiques-et-stress-tests)
10. [Pièges techniques (registre des risques)](#10-pièges-techniques-registre-des-risques)
11. [Stratégie de test et de reproductibilité](#11-stratégie-de-test-et-de-reproductibilité)
12. [Roadmap, portes de validation et chemin critique](#12-roadmap-portes-de-validation-et-chemin-critique)
13. [Plan du rapport final](#13-plan-du-rapport-final)
14. [Références](#14-références)

---

## 0. Thèse et taxonomie du risque

Black-Scholes analytique, l'arbre CRR et Monte Carlo sous GBM sont **un seul modèle et trois méthodes numériques**. Leur convergence vers le même prix est une *porte de correction*, pas une découverte. Le projet étudie trois objets distincts, et le rapport ne les confond jamais :

| Objet | Question | Où on le démontre |
|---|---|---|
| **Risque de méthode numérique** | À quelle vitesse et avec quelle fiabilité chaque méthode converge-t-elle ? Où les implémentations naïves se brisent-elles en silence ? | §3, §4, §5 |
| **Risque de mesure de risque** | Pourquoi deux mesures raisonnables (VaR delta-normale, VaR par réévaluation complète) donnent-elles des réponses radicalement différentes sur la même position ? | §6 |
| **Risque de modèle** | Que se passe-t-il quand les *hypothèses* changent (volatilité stochastique, sauts) alors que le prix des vanilles ATM reste identique ? | §7 |

**Principe directeur :** le risque d'un nombre est une propriété du nombre, pas du code. Un prix Monte Carlo sans erreur standard, un Grec sans analyse biais/variance, une VaR sans intervalle de confiance ni backtest ne sont pas des résultats. L'architecture doit rendre ces omissions difficiles, idéalement impossibles à la compilation.

---

## 1. Ce qui change par rapport à la v1

La v1 est solide sur les mathématiques (toutes ses valeurs de référence ont été reconfirmées). La v2 corrige ou resserre les points suivants :

| # | v1 | v2 | Pourquoi |
|---|---|---|---|
| 1 | `std::mt19937_64` par thread + `std::normal_distribution` | **Philox 4×32-10** (RNG à compteur) + **inverse de la loi normale maison** (Wichura AS241) | `std::normal_distribution` et `std::uniform_real_distribution` sont *définis par l'implémentation* : même graine, résultats différents sous libstdc++, libc++ et MSVC. Incompatible avec une CI 3 compilateurs et un rapport chiffré. |
| 2 | `PriceResult` / `double` | Type **`Estimate`** obligatoire (valeur, erreur standard, biais de discrétisation, n) | Rend impossible la circulation d'un prix MC sans barre d'erreur ; sépare erreur statistique et biais de schéma (critique pour Heston). |
| 3 | Pricer MC avec graine interne | **Pricer stochastique = fonction pure de (marché, graine)** | Les CRN deviennent triviaux et non optionnels. |
| 4 | `std::mdspan`, `std::expected`, `consteval norm_cdf` listés en C++20 | Retirés ou remplacés | `mdspan` et `expected` sont C++23 ; `std::erfc` n'est pas `constexpr` en C++20. Décision actée : **on reste en C++20 strict** (vue 2-D maison, type `Result` maison au lieu de `std::expected`). Pas de bascule C++23 — ce serait de la dérive de périmètre pour un gain marginal. |
| 5 | `-ffast-math` recommandé pour la vectorisation | **Interdit** dans les builds de référence ; autorisé seulement dans un preset `bench-fast` isolé | Casse les tests de NaN, la sommation compensée et la reproductibilité bit-à-bit. |
| 6 | Test MC à 3×SE | **4×SE** avec graine figée | 3σ échoue 0,27 % du temps par test ; avec des dizaines de tests statistiques, les faux positifs deviennent réguliers. |
| 7 | Grecs MC : FD+CRN, pathwise, LR sur un call | **Matrice payoff × estimateur** incluant un **digital** | Le résultat le plus fort est l'échec *silencieux* du pathwise sur un digital (dérivée nulle p.s.) et l'échec des CRN sur payoffs discontinus. |
| 8 | Heston en « stretch » | **Phase à part entière** (§7), mais **explicitement la première chose à couper** si le temps manque (§12.1) | Sans hypothèses alternatives, le mot « modèle » dans *model risk* n'est pas mérité — mais un projet étudiant à temps partiel ne doit pas sacrifier §1–§6 (le cœur, déjà différenciant) pour Heston. |
| 9 | Stress historique = phase thèse | Stress historique **fusionné avec la VaR historique et le backtest** (§6, §9) | Même infrastructure de données ; évite une phase redondante. |
| 10 | Pas de harnais d'expériences | **Phase d'infrastructure dédiée** : chaque figure = un exécutable + CSV + métadonnées | Le rapport doit être régénérable en une commande. |
| 11 | Comparaison des techniques de variance par l'erreur standard | Comparaison par **efficacité = 1 / (variance × temps CPU)** | Une réduction de variance qui coûte plus qu'elle ne rapporte est une régression. |

---

## 2. Architecture C++20

### 2.1 Le problème réel : un produit cartésien partiellement valide

Le moteur manipule **Modèle × Payoff × Méthode**. Toutes les combinaisons n'existent pas : pas d'arbre recombinant pour Heston, pas de formule fermée pour une barrière sous Heston, pas de pathwise utile pour un digital. L'architecture encode cette validité au lieu de la découvrir à l'exécution.

### 2.2 Arborescence

```
RiskEngine-CPP/
├── CMakeLists.txt / CMakePresets.json   # debug, release, asan-ubsan, bench-fast
├── include/riskengine/
│   ├── core/        market.hpp (MarketState), estimate.hpp, units.hpp (typedefs forts),
│   │                rng/philox.hpp, rng/normal_icdf.hpp, stats/welford.hpp
│   ├── models/      gbm.hpp, heston.hpp, merton.hpp        -> concept PathModel
│   ├── payoffs/     vanilla.hpp, digital.hpp, barrier.hpp, asian.hpp
│   ├── methods/     analytic/, lattice/ (crr, leisen_reimer, bbs),
│   │                montecarlo/ (engine, variance_reduction, qmc), fourier/ (heston_cf)
│   ├── greeks/      finite_difference.hpp, pathwise.hpp, likelihood_ratio.hpp, dual.hpp
│   └── risk/        portfolio.hpp, var.hpp, expected_shortfall.hpp, backtest.hpp
├── src/                 # instanciations explicites (extern template)
├── experiments/         # 1 exécutable = 1 figure ou table du rapport
├── tests/               # unitaires, invariants, statistiques
├── bench/               # Google Benchmark
├── data/raw/            # données historiques figées (snapshot versionné)
├── data/results/        # CSV + .meta.json générés
├── tools/               # make_figures.py, check_report_drift.py
└── docs/                # model_risk_report.md, conventions.md, figures/
```

### 2.3 Types fondamentaux

```cpp
// Typedefs forts : tuent la classe de bugs σ vs σ², taux en % vs décimal.
struct Spot     { double v; };
struct Strike   { double v; };
struct Vol      { double v; };   // décimal, par an
struct Rate     { double v; };   // décimal, composition continue
struct Maturity { double v; };   // années (convention documentée dans conventions.md)

struct MarketState {
    Spot spot; Rate rate; Rate div; Vol vol;
};

// Type de retour universel. Aucun prix ne circule en double nu.
struct Estimate {
    double value;
    double std_error      = 0.0;  // erreur statistique (0 si déterministe)
    double discretization = 0.0;  // biais de schéma connu ou estimé
    std::size_t samples   = 0;
};
```

### 2.4 Concepts

```cpp
template <class P>
concept Pricer = requires(const P& p, const MarketState& m) {
    { p.price(m) } -> std::same_as<Estimate>;
};

// Un pricer stochastique est une fonction pure de (marché, clé de graine).
// Aucun état RNG mutable interne : les CRN s'obtiennent en rappelant avec la même clé.
template <class P>
concept StochasticPricer = requires(const P& p, const MarketState& m, SeedKey k) {
    { p.price(m, k) } -> std::same_as<Estimate>;
};

// Dynamique stochastique, découplée des solveurs.
template <class M>
concept PathModel = requires(const M& m, typename M::State s, double dt,
                             std::span<const double> z) {
    { M::factors }         -> std::convertible_to<std::size_t>;
    { m.initial_state() }  -> std::same_as<typename M::State>;
    { m.step(s, dt, z) }   -> std::same_as<typename M::State>;
};

template <class F> concept TerminalPayoff = requires(const F& f, double s) {
    { f(s) } -> std::convertible_to<double>;
};
template <class F> concept PathPayoff = requires(const F& f, std::span<const double> path) {
    { f(path) } -> std::convertible_to<double>;
};

// Validité des combinaisons encodée à la compilation.
template <class Model, class Payoff> concept ClosedForm     = /* traits */;
template <class Model>               concept Recombining    = /* traits */;
template <class Payoff>              concept LipschitzPayoff = /* traits : pathwise admissible */;
```

Chaque type concret est vérifié dans les tests par `static_assert(Pricer<...>)`, afin que l'erreur de concept apparaisse au bon endroit et non au fond d'un `std::visit`.

### 2.5 Pourquoi pas d'héritage, et où garder le `std::variant`

- **Chemin chaud :** templates contraints par concepts. `fd_greeks<Pricer P>` est écrit une seule fois, inliné, sans appel virtuel.
- **Couche de comparaison :** `std::variant` + `std::visit`, conservé de la v1 (ensemble fermé, sémantique de valeur, pas d'allocation, exhaustivité vérifiée).
- **Limite de la v1 :** dès que Heston et Merton entrent, un variant de *méthodes* ne suffit plus. Le variant porte sur des **scénarios pré-validés** `Scenario<Model, Payoff, Method>`, construits par une fabrique qui refuse les combinaisons invalides à la compilation.
- **Temps de compilation :** instanciations explicites dans `src/`, `extern template` dans les headers.
- **Piège à anticiper (nouveau en v2) :** un `static_assert` de concept mal formé au fond d'une fabrique template peut produire un mur d'erreurs illisible. Chaque contrainte de fabrique doit porter un message explicite (`static_assert(ClosedForm<Model,Payoff>, "pas de formule fermée pour cette combinaison")`), testé isolément avant d'être composé.

### 2.6 Différentiation automatique forward par templates

Les modèles sont templatés sur le scalaire : `GBM<Real>`. Instancié avec `Dual<double>` (≈ 80 lignes), le moteur MC produit **les Grecs pathwise automatiquement**, ce qui fournit une vérification croisée indépendante des dérivées pathwise écrites à la main. L'AAD (mode adjoint) reste hors chemin critique — à ne considérer qu'après la phase 8, jamais avant.

### 2.7 Ce que l'on prend et ce que l'on laisse de QuantLib

- **Pris :** la séparation `Instrument` / `PricingEngine` — exactement l'exigence « un instrument, plusieurs moteurs ».
- **Laissé :** `shared_ptr` omniprésents, `Handle<>`, hiérarchies virtuelles profondes dans les boucles chaudes. Savoir expliquer *pourquoi* on fait différemment en C++20 est plus fort que l'une ou l'autre approche seule.

---

## 3. Fondations mathématiques et valeurs de référence

### 3.1 Black-Scholes avec dividende continu

```
d1 = [ln(S/K) + (r − q + σ²/2)T] / (σ√T),   d2 = d1 − σ√T
C  = S e^{−qT} N(d1) − K e^{−rT} N(d2)
P  = K e^{−rT} N(−d2) − S e^{−qT} N(−d1)
```

**Vecteur canonique ✅** (S = K = 100, r = 5 %, q = 0, σ = 20 %, T = 1) :

| Quantité | Valeur |
|---|---|
| d1, d2 | 0,35 ; 0,15 |
| Call | 10,450584 |
| Put | 5,573526 |
| Delta (call) | 0,636831 |
| Gamma | 0,018762 |
| Vega (par 1,00 de vol) | 37,524035 → **0,375 par point de vol** |
| Theta (call, par an) | −6,414028 → **−0,01757 par jour calendaire** |
| Rho (par 1,00 de taux) | 53,232482 |

**Invariants testés** (au-delà de « prix ≈ attendu ») :

- parité call-put à ~1e-14 ;
- bornes : max(0, S e^{−qT} − K e^{−rT}) ≤ C ≤ S e^{−qT} ;
- monotonie : ∂C/∂K < 0, ∂C/∂σ > 0, ∂C/∂S > 0 ;
- convexité en strike (butterfly ≥ 0) : absence d'arbitrage statique ;
- queues : N(x) calculé par `0.5 * std::erfc(-x / std::sqrt(2.0))`, jamais par `1 − N(−x)`.

**Volatilité implicite :** Jäckel, *Let's Be Rational* (2015) comme cible ; Brent sur un intervalle encadrant en repli ; **jamais Newton seul** (divergence en ailes profondes où vega → 0). Porte : aller-retour prix → vol → prix sur une grille couvrant moneyness 0,5–2 et T de 1 jour à 5 ans.

### 3.2 Arbres binomiaux

```
Δt = T/n,  u = e^{σ√Δt},  d = 1/u,  p = (e^{(r−q)Δt} − d)/(u − d)
```

- **Validité :** p ∈ (0,1) ⇔ Δt < σ²/(r − q)². Assertion obligatoire.
- **Mémoire :** induction arrière en place sur un seul vecteur, O(n) (v1 : 400 Mo → 0,08 Mo à n = 10 000).
- **Convergence non monotone :** l'erreur de CRR change de signe à chaque n (position du strike entre les nœuds terminaux, Leisen-Reimer 1996). Échantillonner seulement des n pairs produit une pente de −1 artificiellement propre. Le rapport trace des **n consécutifs**.
- **Remèdes étudiés :** moyenne (Vₙ + Vₙ₊₁)/2 ; **Leisen-Reimer** (inversion de Peizer-Pratt, n impair, convergence O(1/n²)) ; **BBS** (Broadie-Detemple : BS à l'avant-dernier pas) + extrapolation de Richardson.
- **Grecs gratuits :** delta et gamma lus sur les nœuds des pas 1 et 2 ; theta via un arbre démarrant deux pas avant t = 0.
- **Invariants américains ✅ (n = 20 000) :** put américain 6,090333 > put européen 5,573426 (prime d'exercice anticipé 0,516907) ; call américain sans dividende = call européen exactement (théorème de Merton).

### 3.3 Monte Carlo sous GBM

- Payoff européen : **tirage direct de S_T** en log, sans pas de temps (exact, sans biais de discrétisation).
- Payoff dépendant du chemin : pas de temps en log (exact aux dates d'observation pour GBM) ; le biais de *monitoring discret* d'une barrière est un effet distinct, documenté et corrigible (Broadie-Glasserman-Kou).
- Convergence en O(N^{−1/2}) vérifiée par pente log-log **et** par **couverture empirique** : 1 000 pricings indépendants, ~95 % des IC à 95 % doivent contenir le prix BS. Presque aucun projet étudiant ne fait ce test.

---

## 4. Monte Carlo : échantillonnage et réduction de variance

### 4.1 Métrique

**Efficacité = 1 / (Var × temps CPU)** (Glasserman). Chaque table du rapport donne variance, temps et efficacité relative. La v1 a mesuré, à tirages égaux : antithétique 4×, contrôle S_T 6,9×, antithétique + contrôle **58×** (en équivalent de chemins). Message à conserver : **la réduction de variance bat le parallélisme** (58× par les mathématiques contre ~6× par 8 threads) — à reconfirmer en efficacité CPU.

### 4.2 Techniques, avec leurs échecs

| Technique | Fonctionne quand | Échoue quand (à montrer) |
|---|---|---|
| Antithétique | Payoff monotone en Z (call, put) | Straddle ATM : payoff quasi pair en Z, gain nul voire négatif en efficacité |
| Contrôle S_T, E = S₀e^{(r−q)T} | Calls ITM, forte corrélation | Calls très OTM : corrélation qui s'effondre |
| Contrôle asiatique géométrique | Asiatique arithmétique (formule fermée pour le géométrique, corrélation > 0,99) | — démonstration « spectaculaire » |
| Contrôle BS sous Heston | Même Z, gain modéré | Se dégrade avec le vol-of-vol : mesure indirecte de l'écart de modèle |
| QMC Sobol (Joe-Kuo) + pont brownien + brouillage d'Owen | Faible dimension effective, payoffs lisses | Payoffs discontinus ; sans randomisation, plus d'erreur standard |
| Importance sampling (décalage de moyenne) | Digitals très OTM, queues de VaR | Mauvais choix du décalage : variance qui explose |

**Détails qui comptent :**

- Le β* optimal du contrôle est estimé sur un **pilote indépendant** ; l'estimer sur les mêmes chemins introduit un biais en O(1/N), petit mais à mentionner.
- En QMC, **pas de Box-Muller** : il détruit la structure de basse discrépance. Inverse de la fonction de répartition uniquement.

### 4.3 Générateur et parallélisme reproductible

- **Philox 4×32-10** : sous-flux indexé par (graine, identifiant de bloc), sans état partagé, sans saut d'état.
- **Décomposition en K blocs fixes**, indépendante du nombre de threads. Chaque bloc produit un accumulateur de Welford ; les partiels sont réduits **dans l'ordre des blocs** (fusion de Chan). Résultat **bit-à-bit identique** avec 1 ou 64 threads.
- Porte de la phase 2 : même résultat au bit près avec 1, 2, 8 threads, sur GCC et Clang.
- **Contingence explicite (nouveau en v2) :** Philox + AS241 sont de l'infrastructure neuve, pas des formules connues copiées d'un manuel — c'est la partie du plan la plus susceptible de déraper. Si la phase 2 dépasse 1,5 semaine, replier sur `std::mt19937_64` + inverse normale maison (garde la reproductibilité *par compilateur*, perd la garantie *inter*-compilateur) plutôt que de laisser le retard se propager sur tout le chemin critique. Documenter le repli comme limite connue dans le rapport si utilisé.

---

## 5. Estimation des Grecs sous bruit (section phare)

### 5.1 Protocole

Quatre estimateurs × deux payoffs (**call** : continu avec un pli ; **digital** : discontinu) × trois régimes (ATM, OTM profond, T court). Pour chaque cellule : biais, variance, RMSE, coût, efficacité, en fonction de h et de N.

### 5.2 Résultats théoriques à confirmer empiriquement

| Estimateur | Call (lipschitzien) | Digital (discontinu) |
|---|---|---|
| **FD centrée, graines indépendantes** | Biais O(h²), Var O(1/(N h²)) → h* ∝ N^{−1/6}, RMSE ∝ N^{−1/3}. Explose quand h → 0. | Idem, pire |
| **FD centrée + CRN** | Var O(1/N), indépendante de h → RMSE ∝ N^{−1/2} | Var O(1/(N h)) → h* ∝ N^{−1/5}, RMSE ∝ N^{−2/5}. **Les CRN ne sauvent pas tout.** |
| **Pathwise** | Sans biais, faible variance, une seule passe | **Vaut 0 presque sûrement : converge avec assurance vers une réponse fausse** |
| **Likelihood ratio** | Sans biais, variance plus élevée | Sans biais ; variance ∝ 1/T quand T → 0 |

**Illustration mesurée en v1 (200 000 chemins, delta vrai 0,636831) :** avec h = 0,01, FD sans CRN donne 3,387 (erreur +275 %) ; avec CRN, 0,6345, stable sur quatre ordres de grandeur de h.

**Gamma par FD d'un call, même avec CRN :** le pli du payoff fait que la seconde différence ne contribue que dans une bande de largeur h autour du strike, avec une amplitude en 1/h ; la variance se comporte en O(1/(N h)), exactement comme un delta de digital.

### 5.3 Formules (GBM, Z ~ N(0,1) pilotant S_T)

```
Pathwise   Δ = e^{−rT} E[ 1{S_T > K} · S_T / S₀ ]
           ν = e^{−rT} E[ 1{S_T > K} · S_T (√T·Z − σT) ]
LR         Δ = e^{−rT} E[ payoff · Z / (S₀ σ √T) ]
           Γ = e^{−rT} E[ payoff · (Z² − 1 − Zσ√T) / (S₀² σ² T) ]
Mixte      Γ = LR appliqué à l'estimateur pathwise du delta (variance plus faible que LR pur)
```

v1 : delta pathwise 0,634511 et gamma LR 0,018554 (vrai 0,018762), à reproduire avec erreurs standard.

### 5.4 Arrondi flottant sur les pricers déterministes

Sur BS analytique, tracer l'erreur de la FD en fonction de h : **courbe en V** (troncature O(h²) à droite, arrondi O(ε/h) pour le delta et O(ε/h²) pour le gamma à gauche). h* ≈ S·ε^{1/3} pour le delta, ≈ S·ε^{1/4} pour le gamma. Bump **relatif**, jamais absolu.

### 5.5 Livrable

Une **matrice de recommandations** payoff × régime × estimateur, et une règle de décision opérationnelle pour une équipe de validation.

---

## 6. Mesures de risque sur positions non linéaires

### 6.1 La démonstration centrale

Position : **short straddle ATM, delta-hedgé**. Delta ≈ 0, donc **VaR delta-normale ≈ 0**, alors que la perte réelle ½·Γ·ΔS² est grande et asymétrique.

**Aperçu ✅** (straddle 30 jours, σ = 20 %, r = 5 %, horizon 1 jour, 10⁶ scénarios, par straddle) : VaR delta-normale = **0** ; VaR 99 % par réévaluation complète ≈ **0,64** ; ES 97,5 % ≈ **0,66**. À régénérer et à étendre en phase 6.

### 6.2 Progression des méthodes

1. **Delta-normale** : aveugle au gamma par construction.
2. **Delta-gamma normale** (moments d'une forme quadratique projetés sur une normale) : amélioration partielle ; la vraie loi est de type χ² décentrée, asymétrique, et la queue reste sous-estimée. **Cornish-Fisher** améliore sans résoudre.
3. **Réévaluation complète** Monte Carlo et historique : capture l'optionalité.
4. **Chocs de volatilité conjoints** : le short straddle est aussi short vega ; une VaR sur le seul spot rate ce risque.

### 6.3 Conventions (source de la majorité des bugs)

- VaR rapportée comme **perte positive** ; 99 % = 1er percentile du P&L.
- Estimateur de quantile : index ⌈αN⌉ sur les pertes triées ; `std::nth_element` (O(N)).
- Mise à l'échelle √10 : suppose des rendements i.i.d., **fausse** sous clustering de volatilité ; à démontrer, pas à utiliser en silence.
- **IC bootstrap sur toute VaR** : la variance de l'estimateur de quantile ∝ α(1−α)/(N f(q)²), énorme en queue (v1 : 5,6 % d'erreur à N = 1 000).

### 6.4 Expected Shortfall : non optionnel

- **Non-sous-additivité de la VaR ✅** : deux obligations indépendantes, défaut 4 %, perte 100 : VaR95(A) = VaR95(B) = 0 mais VaR95(A+B) = 100. La diversification *augmente* la VaR (Artzner et al., 1999).
- **FRTB** (BCBS d457, 2019) : ES 97,5 % remplace la VaR 99 %. Sous normalité, ES 97,5 % = 2,3378 ≈ VaR 99 % = 2,3263.

| α | VaR (N(0,1)) ✅ | ES ✅ |
|---|---|---|
| 95 % | 1,6449 | 2,0627 |
| 97,5 % | 1,9600 | 2,3378 |
| 99 % | 2,3263 | 2,6652 |

### 6.5 Backtesting

- **Kupiec POF** (couverture inconditionnelle) et **Christoffersen** (indépendance des dépassements).
- **Feux de Bâle** sur 250 jours : vert 0–4 dépassements, jaune 5–9, rouge ≥ 10.
- Résultat attendu : sur données réelles, les dépassements sont groupés en crise ; c'est la signature des queues épaisses et du clustering que GBM ignore.

---

## 7. Risque de modèle au-delà de GBM

**Cadrage (nouveau en v2) : cette phase est explicitement la première extension à couper si le temps manque (§12.1).** Le cœur du projet (§1–§6) est déjà différenciant et complet sans elle. Ne pas commencer la phase 7 tant que la phase 6 n'est pas validée et rédigée.

### 7.1 Modèles

- **Heston** : simulation par schéma **QE d'Andersen (2008)** ; Euler à troncature complète comme contre-exemple ; prix semi-analytique par fonction caractéristique dans la formulation d'**Albrecher et al. (2007)** (évite la coupure du logarithme complexe). Condition de Feller 2κθ > ξ² : étudier les deux régimes.
- **Merton (sauts log-normaux)** : prix fermé comme série pondérée de prix BS ; simulation exacte (Poisson + normale).

### 7.2 Trois expériences

1. **Même prix ATM, prix d'exotiques divergents.** Calibrer BS, Heston et Merton sur le même prix ATM (ou le même petit smile), puis pricer un call OTM, un digital et une barrière up-and-out. L'écart entre modèles *est* la mesure du risque de modèle. Les barrières sont notoirement sensibles au vol-of-vol même quand les vanilles coïncident.
2. **Couverture dans le mauvais monde.** Delta-hedge avec le delta BS alors que le monde simulé est Heston ou Merton ; tracer la distribution du P&L de couverture en fonction de la fréquence de rebalancement. Comparer à P&L ≈ ½ Γ S² (σ²_réalisée − σ²_implicite) dt. Sous Merton, un P&L résiduel non couvrable subsiste quelle que soit la fréquence.
3. **Décomposition des erreurs sous Heston.** Biais de discrétisation (↓ avec dt) contre erreur statistique (↓ avec N) ; montrer que QE élimine l'essentiel du premier. C'est ici que le champ `Estimate::discretization` prend tout son sens.

---

## 8. Ingénierie de performance

Section volontairement secondaire dans le rapport ; utile, mais ce n'est pas la thèse.

- **Threads :** `std::jthread` + pool de blocs (ou OpenMP). Éviter `std::async` (politique de lancement définie par l'implémentation) et `std::execution::par` (support inégal : libstdc++ dépend de TBB, libc++ incomplet).
- **False sharing :** accumulateurs locaux au thread, écrits une fois ; `alignas(std::hardware_destructive_interference_size)` quand un tableau partagé est inévitable. v1 mesurait un ralentissement de ~9,7× ; à re-mesurer et à montrer comme étude de profilage.
- **SoA et lots de normales** (1 024 à 4 096) pour permettre la vectorisation ; `exp` domine la boucle interne.
- **Arbres :** puissances de u précalculées au lieu de `std::pow` dans la boucle interne.
- **Cibles indicatives (v1, à reconfirmer) :** BS < 100 ns ; arbre n = 1 000 < 5 ms ; MC 1M chemins mono-thread < 20 ms ; efficacité parallèle à 8 threads > 70 %.
- **Rigueur de mesure :** Google Benchmark, médiane et dispersion, modèle de CPU, compilateur et drapeaux consignés. Expliquer l'écart à la linéarité (hyperthreading, baisse de fréquence en charge, bande passante mémoire) vaut mieux que prétendre à un speedup linéaire.

---

## 9. Données historiques et stress tests

### 9.1 Sources (à revérifier avant la phase 6)

| Source | Usage | Remarques |
|---|---|---|
| **FRED** | Taux sans risque historiques (DGS3MO, DGS1), VIXCLS | Officiel ; clé API gratuite. Source privilégiée pour r et la vol implicite proxy. |
| **Stooq** | Historique OHLCV d'indices et d'ETF, export CSV | Sans inscription ; bon repli. |
| **Yahoo Finance** (`yfinance`) | Historique OHLCV | API non officielle, casse parfois. |
| **EODHD** | Historique long d'indices | Quotas du niveau gratuit à vérifier. |

**Règle de reproductibilité :** les données sont téléchargées une fois, **figées dans `data/raw/`** avec leur date d'extraction et leur somme de contrôle. Aucune expérience ne télécharge à l'exécution.

**Limite assumée :** les chaînes d'options historiques (2008, 2020) ne sont pas disponibles gratuitement. Conception honnête : prix réels du sous-jacent + VIX comme proxy de la vol implicite, pour pricer une option *hypothétique* ; limite énoncée explicitement dans le rapport.

### 9.2 Scénarios

Octobre 1987, septembre–octobre 2008, février 2018 (« Volmageddon »), mars 2020. Chocs **conjoints** spot + vol + taux, jamais un facteur à la fois. Réutilise l'infrastructure de la VaR historique et du backtest (§6).

---

## 10. Pièges techniques (registre des risques)

| Piège | Symptôme | Parade |
|---|---|---|
| Distributions `std::` définies par l'implémentation | Chiffres différents selon le compilateur | Philox + inverse de la normale maison |
| Réduction parallèle non déterministe | Résultat qui dépend du nombre de threads | Blocs fixes, réduction ordonnée |
| `-ffast-math`, contraction FMA | NaN non détectés, dérive entre compilateurs | Interdit en référence ; tolérances relatives dans les tests |
| Annulation catastrophique en FD | Gamma bruité ou faux pour h petit | Bump relatif, h* documenté, courbe en V dans le rapport |
| Queue de N(x) via 1 − N(−x) | Précision nulle en queue | `erfc` |
| T → 0, σ → 0 | Division par zéro, NaN | Branches explicites vers la valeur intrinsèque, tests dédiés |
| p hors de (0,1) dans CRR | Prix plausibles mais faux, arbitrage | Assertion |
| Variance négative sous Heston-Euler | NaN ou biais massif | Troncature complète ou QE |
| Coupure du logarithme complexe (Heston CF) | Prix semi-analytiques faux pour T long | Formulation d'Albrecher et al. |
| Box-Muller en QMC | Perte de la basse discrépance | Inverse de la fonction de répartition |
| Conventions (252/365, vega par 1 % ou 1,00, theta par jour ou an) | Erreurs d'un facteur 100 ou 365 | `docs/conventions.md` dès la phase 1, typedefs forts |
| Tests MC à tolérance fixe | Tests instables ou trop laxistes | |MC − ref| < 4·SE, graine figée |
| Erreurs de concept illisibles | Diagnostic au fond d'un `std::visit` | `static_assert` par type dans les tests |
| Explosion du temps de compilation | Builds lents | Instanciations explicites, `extern template` |
| Dérive de périmètre | Projet jamais fini | **Aucune fonctionnalité sans figure prévue dans le rapport** |

---

## 11. Stratégie de test et de reproductibilité

**Trois niveaux de tests :**

1. **Unitaires et valeurs de référence** : vecteurs canoniques (§3.1), valeurs publiées pour l'américain et Heston.
2. **Invariants** : parité, bornes, monotonie, convexité, américain ≥ européen, call américain sans dividende = européen.
3. **Statistiques** : |MC − ref| < 4·SE ; pentes de convergence dans un intervalle ; couverture des IC ; moments et Kolmogorov-Smirnov du générateur.

**Harnais d'expériences :** chaque exécutable de `experiments/` écrit `data/results/<id>.csv` et `<id>.meta.json` (SHA git, compilateur, drapeaux, graine, N, date). `tools/make_figures.py` produit `docs/figures/`. `make report` régénère tout.

**Le rapport comme test :** un job CI hebdomadaire relance les expériences et échoue si un chiffre cité dans le rapport dérive au-delà de sa tolérance statistique.

**CI :** GCC, Clang, MSVC ; `-Wall -Wextra -Werror` ; job ASan/UBSan ; clang-format et clang-tidy.

---

## 12. Roadmap, portes de validation et chemin critique

Chaque phase a une **porte de validation** et **alimente une section du rapport**. Une phase n'est terminée que lorsque sa section est rédigée.

| Phase | Contenu | Porte de validation | Rapport | Durée indicative (temps partiel) |
|---|---|---|---|---|
| **0** ✅ | Scaffold : CMake, CI, concept `Pricer` | CI verte sur 3 compilateurs | — | fait |
| **1** | BS avec dividende, Grecs fermés, vol implicite, `conventions.md` | Vecteur canonique à 1e-9, parité à 1e-14, aller-retour vol implicite sur toute la grille, cas limites sans NaN | §3 | 1 sem. |
| **2** | Philox, inverse normale, `Estimate`, Welford, réduction ordonnée, harnais d'expériences | Tests statistiques du RNG ; résultat bit-à-bit identique à 1/2/8 threads et sur 2 compilateurs | §2 | 1 sem. (repli à 1,5 sem. max — voir §4.3) |
| **3** | CRR, moyenne, Leisen-Reimer, BBS + Richardson, américain, Grecs sur nœuds | Pentes −1 (CRR) et −2 (LR) ; oscillation visible sur n consécutifs ; invariants américains | §4 | 1 sem. |
| **4** | Moteur MC générique `PathModel`, antithétique, contrôles, QMC randomisé, asiatique | |MC − BS| < 4·SE ; pente −0,5 ; couverture des IC ≈ 95 % ; table d'efficacité | §5 | 1,5 sem. |
| **5** | Grecs : FD naïve, FD+CRN, pathwise, LR, mixte, vérification par `Dual` | Matrice payoff × estimateur complète ; FD+CRN stable sur le call, instable sur le digital ; pathwise digital = 0 démontré | §6 | 1,5 sem. |
| **6** | VaR/ES (delta-normale, delta-gamma, Cornish-Fisher, MC, historique), bootstrap, backtests, stress historiques | Table normale analytique retrouvée ; contre-exemple de sous-additivité ; VaR delta-normale ≈ 0 sur le short straddle contre VaR réelle significative | §7 | 2 sem. |
| **7** *(optionnelle — voir §12.1)* | Heston (QE + CF), Merton, trois expériences de risque de modèle | Prix Heston CF contre MC-QE dans 4·SE ; Merton série contre MC ; distributions de P&L de couverture | §8 | 2 sem. |
| **8** | Benchmarks, rédaction finale, reproduction intégrale | `make report` régénère toutes les figures ; relecture par un non-spécialiste | §9–11 | 1 sem. |

**Chemin critique : 1 → 2 → 4 → 5 → 6 → 8.** (La phase 7 n'est *pas* sur le chemin critique — voir §12.1.)

- La phase 3 est **hors chemin critique** et parallélisable avec la 4 (sauf si l'on veut l'américain sous Heston en phase 7).
- **Le risque de planning principal est la phase 2 bâclée** : un générateur non reproductible invalide rétroactivement tous les chiffres du rapport. Voir le plan de repli en §4.3.
- Total indicatif : **~9 semaines à temps partiel pour un projet complet et rédigé (phases 0–6 + 8), ~11–12 semaines si la phase 7 est incluse.**

### 12.1 Descoping — quoi couper si le temps manque

Dans l'ordre, la première chose à sacrifier est toujours la plus récemment ajoutée au périmètre, jamais le cœur validé :

1. **Phase 7 (Heston/Merton) en entier.** Le projet reste complet et cohérent sans elle : §1–§6 démontrent déjà les trois volets de la thèse (méthode numérique, mesure de risque) sauf le risque de modèle au sens strict. Le README doit alors dire explicitement *"risque de méthode numérique et de mesure de risque démontrés ; risque de modèle au sens strict identifié comme extension future, non traité"* — honnête et toujours un signal fort.
2. **Extensions du §12 (tableau) :** Longstaff-Schwartz, AAD adjoint, SIMD explicite — jamais commencées avant la phase 8.
3. **Dans la phase 6 :** garder VaR/ES + sous-additivité + le short straddle (le résultat phare) ; le backtesting complet (Kupiec + Christoffersen + feux de Bâle) peut se réduire à Kupiec seul si nécessaire.
4. **Ne jamais couper :** phase 1 (fondation), phase 2 (reproductibilité — tout en dépend), la matrice de Grecs de la phase 5 (résultat le plus différenciant du projet).

**Extensions, classées par rapport valeur/effort (à ne considérer qu'après la phase 8) :**

| Extension | Valeur | Effort |
|---|---|---|
| Correction de monitoring discret des barrières (BGK) | Élevée (lien §4/§7) | Faible |
| Liaisons Python (nanobind) pour les expériences | Moyenne-élevée | Faible |
| Longstaff-Schwartz (américain par MC) | Élevée | Élevé |
| AAD mode adjoint | Élevée | Élevé |
| Volatilité locale (Dupire) comme 4ᵉ hypothèse | Très élevée pour la thèse | Élevé |
| SIMD explicite | Moyenne, seulement si mesurée | Élevé |

---

## 13. Plan du rapport final

`docs/model_risk_report.md`, rédigé comme un papier de validation de modèles de niveau institutionnel.

**Résumé exécutif** — trois ou quatre découvertes chiffrées, en un paragraphe, avec leurs intervalles de confiance.

1. **Introduction et taxonomie du risque**
   1.1 Risque de méthode numérique, risque de mesure, risque de modèle
   1.2 Périmètre : ce que le rapport démontre et ce qu'il ne prétend pas démontrer
2. **Cadre, conventions et reproductibilité**
   2.1 Dynamiques, paramètres de référence, conventions d'unités
   2.2 Protocole : générateur, graines, métrique d'efficacité, régénération
3. **Vérité terrain : Black-Scholes analytique**
   3.1 Validation, invariants, volatilité implicite
   3.2 Arrondi flottant et courbe en V des différences finies
4. **Arbres : convergence et pathologies**
   4.1 Oscillations pair/impair de CRR et leur origine
   4.2 Leisen-Reimer, BBS, Richardson : ordres mesurés
   4.3 Exercice anticipé
5. **Monte Carlo : convergence et réduction de variance**
   5.1 Taux N^{−1/2} et couverture des intervalles
   5.2 Table d'efficacité, échecs compris
   5.3 Quasi-Monte Carlo randomisé
6. **Estimation des Grecs sous bruit** *(section phare)*
   6.1 Pourquoi la FD naïve échoue : analyse biais/variance
   6.2 CRN : ce qu'ils corrigent et ce qu'ils ne corrigent pas
   6.3 Pathwise, LR, estimateur mixte ; l'échec silencieux sur le digital
   6.4 Matrice de recommandations
7. **Mesures de risque sur positions non linéaires**
   7.1 La cécité de la VaR delta-normale sur un short straddle
   7.2 Delta-gamma et Cornish-Fisher
   7.3 Réévaluation complète, historique, Expected Shortfall, sous-additivité
   7.4 Backtesting et scénarios de stress historiques
8. **Risque de modèle au-delà de GBM** *(si la phase 7 a été réalisée — sinon, section courte expliquant le périmètre non traité, voir §12.1)*
   8.1 Calibration commune, exotiques divergents
   8.2 P&L de couverture dans un monde mal spécifié
   8.3 Biais de discrétisation contre erreur statistique
9. **Notes d'ingénierie et performance** *(courte)*
10. **Limites et travaux futurs**
11. **Conclusion : recommandations opérationnelles** pour un desk ou une équipe de validation

**Annexes :** A. dérivations des estimateurs pathwise et LR ; B. tables de paramètres ; C. valeurs de référence et sources ; D. procédure de reproduction.

**Règles éditoriales :**

- Chaque figure a une légende autonome : *ce qu'elle montre*, puis *pourquoi c'est important*.
- Chaque chiffre renvoie au CSV qui l'a produit.
- Chaque estimation stochastique est donnée avec son erreur standard ou son IC.
- Chaque affirmation négative (« X échoue ») est accompagnée du mécanisme et de la parade.

---

## 14. Références

- Hull, *Options, Futures, and Other Derivatives* (numéros de chapitres selon l'édition).
- Glasserman, *Monte Carlo Methods in Financial Engineering* (2003) — ch. 4 (réduction de variance), ch. 7 (Grecs).
- Joshi, *C++ Design Patterns and Derivatives Pricing*.
- Leisen & Reimer (1996), *Binomial models for option valuation — examining and improving convergence*, Applied Mathematical Finance 3(4).
- Broadie & Detemple (1996), *American option valuation: new bounds, approximations, and a comparison of existing methods*, Review of Financial Studies 9(4).
- Jäckel (2015), *Let's Be Rational*, Wilmott.
- Salmon, Moraes, Dror, Shaw (2011), *Parallel random numbers: as easy as 1, 2, 3* (Philox / Random123).
- Joe & Kuo (2008), *Constructing Sobol sequences with better two-dimensional projections*.
- Heston (1993), *A closed-form solution for options with stochastic volatility*, RFS 6(2).
- Albrecher, Mayer, Schoutens, Tistaert (2007), *The little Heston trap*, Wilmott.
- Andersen (2008), *Simple and efficient simulation of the Heston stochastic volatility model*, Journal of Computational Finance 11(3).
- Merton (1976), *Option pricing when underlying stock returns are discontinuous*, JFE 3.
- Artzner, Delbaen, Eber, Heath (1999), *Coherent Measures of Risk*, Mathematical Finance 9(3).
- Kupiec (1995), POF test ; Christoffersen (1998), conditional coverage.
- BCBS d457 (2019), *Minimum capital requirements for market risk* (FRTB).
- Longstaff & Schwartz (2001), *Valuing American Options by Simulation*, RFS 14(1).
- Code source de QuantLib (`ql/instrument.hpp`, `ql/pricingengine.hpp`) — pour l'abstraction, pas pour le style.
