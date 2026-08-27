/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSeriesNatural
import ABC3.Found.GaloisRep.TateOrigin

/-!
# Tate 一意化の**座標は自然**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★段 2c-2 の到達点

`galois-equivariant-tate-uniformization` の段 2c は
「`tatePtPair` の自然性（級数と `σ` の交換）」であり、同変性の**本体**である。

★本ファイルで **`tatePtPair` の座標そのもの**が自然であることが取れる:

> **`σ_K (X_K(a,w,q)) = X_K(σa, σw, σq)`**、`Y` も同様

★★`tatePtPair a w q … = Point.some (tateXK a w q hq) (tateYK a w q hq) …` なので、
残るのは `Point.some` の合同（段 2c-3）だけである。

## ★★★★★★積み上げの鎖

| 段 | 補題 | ファイル |
|---|---|---|
| 2c-1 | `evalAdic_map`（冪級数の adic 値） | `AdicEvalNatural.lean` |
| 2c-2 前半 | `adicSum_map`・`tateXtail_map`・`tateYtail_map` | `TateSeriesNatural.lean` |
| 2c-2 後半 | `tateXpairE_map`・**`tateXK_map`** | ★本ファイル |
| 2c-3 | `tatePtPair` の自然性 | 未 |

★★★どの段も「**一意 ⟹ 自然**」か「環演算」だけである。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★`σ` は `σR : R →+* R` と `σK : K →+* K` の**対**として受け、
両立 `σK ∘ algebraMap = algebraMap ∘ σR` と `σ I ⊆ I` を仮説に置く。
★★体の自己同型から来る場合、`σ I ⊆ I` は付値の拡大の一意性から出るが、
本ファイルはそれを仮説にしている（段 3 の担当）。
-/

namespace ABC3.Found.GaloisRep

/-! ## ★★★★分母を払った座標の自然性 -/

/-- ★★★★★★★**`(1−a)²·X` は自然**。 -/
theorem tateXpairE_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (a w q : R)
    (hq : q ∈ I) (hq' : σ q ∈ I) (hw : w ∈ I) :
    σ (tateXpairE a w q hq) = tateXpairE (σ a) (σ w) (σ q) hq' := by
  unfold tateXpairE
  rw [map_add, map_mul, map_pow, map_sub, map_one, map_sub, map_add, map_add,
    map_mul, map_ofNat,
    tateXtail_map σ hσI a q hq hq', tateXtail_map σ hσI w q hq hq',
    tateXterm_map σ (I := I) hw,
    evalAdic_map σ hσI _ q hq hq']

/-- ★★★★★★★**`(1−a)³·Y` は自然**。 -/
theorem tateYpairE_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (a w q : R)
    (hq : q ∈ I) (hq' : σ q ∈ I) (hw : w ∈ I) :
    σ (tateYpairE a w q hq) = tateYpairE (σ a) (σ w) (σ q) hq' := by
  unfold tateYpairE
  simp only [map_add, map_sub, map_mul, map_pow, map_one]
  rw [tateYtail_map σ hσI a q hq hq', tateXtail_map σ hσI w q hq hq',
    tateYtail_map σ hσI w q hq hq',
    tateXterm_map σ (I := I) hw, tateYterm_map σ (I := I) hw,
    evalAdic_map σ hσI _ q hq hq']

/-! ## ★★★★★★★★★★`K` の水準の座標 -/

/-- ★★★★★★★★★★**Tate 一意化の `X` 座標は自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> **`σ_K (X_K(a,w,q)) = X_K(σa, σw, σq)`**

★★`tateXK = algebraMap R K (tateXpairE …) · (algebraMap R K (1−a))⁻¹²` なので、
`R` の中の自然性（上）と `σ_K` が体準同型であること（`map_inv₀`）だけで出る。 -/
theorem tateXK_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hq' : σ q ∈ I) (hw : w ∈ I) :
    σK (tateXK a w q hq : K) = tateXK (σ a) (σ w) (σ q) hq' := by
  unfold tateXK
  rw [map_mul, map_pow, map_inv₀, hcompat, hcompat, map_sub, map_one,
    tateXpairE_map σ hσI a w q hq hq' hw]

/-- ★★★★★★★★★★**Tate 一意化の `Y` 座標は自然**。 -/
theorem tateYK_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σ r))
    (a w q : R) (hq : q ∈ I) (hq' : σ q ∈ I) (hw : w ∈ I) :
    σK (tateYK a w q hq : K) = tateYK (σ a) (σ w) (σ q) hq' := by
  unfold tateYK
  rw [map_mul, map_pow, map_inv₀, hcompat, hcompat, map_sub, map_one,
    tateYpairE_map σ hσI a w q hq hq' hw]

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**段 2c-2 の到達点**であって、
`tatePtPair`（`Point.some` の合同、段 2c-3）ではない。 -/

def tateXpairE_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分母を払った X の自然性)",
    sectionId := "genell-def-3-3" }

def tateXK_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X 座標の自然性)",
    sectionId := "genell-def-3-3" }

def tateYK_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Y 座標の自然性)",
    sectionId := "genell-def-3-3" }

def tateXK_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tateXtail_map / evalAdic_map(級数の自然性——段 2c-1・2c-2 前半)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXtail_map") 15,
    .citation "[mathlib]" "map_inv₀(体準同型は逆元と可換)"
      (.inMathlib "map_inv₀") 15,
    .implicitStep
      ("★★★★これで**tatePtPair の座標そのもの**が自然であることが取れた。" ++
       "★残るのは Point.some の合同(段 2c-3)だけである") 15,
    .implicitStep
      ("★逸脱: σ を σR : R →+* R と σK : K →+* K の**対**として受け、" ++
       "両立と σ I ⊆ I を仮説に置いた。" ++
       "★★体の自己同型から来る場合 σ I ⊆ I は付値の拡大の一意性から出る(段 3)") 15 ]

end ABC3.Found.GaloisRep
