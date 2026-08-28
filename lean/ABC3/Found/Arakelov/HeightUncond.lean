/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegAPicM
import ABC3.Found.Arakelov.PullbackPic
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★**計量表示の高さ（無条件）** `ht_M̄(x_F) = deg_F(x_F^* M̄)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

`MetricHeight.lean` の `htMetric` は `SzpData`（[Szp] Prop 1.1 の同型）を
**仮定として担いでいた**。
★`§9-782` で `deg_F : APic(Spec 𝓞_F) → ℝ` が**無条件に**作れたので、
本ファイルはその仮定を外した版を置く。

    `ht_M̄(x_F) ≝ deg_F(x_F^* M̄)`,  `ht_{M̄ ⊗ N̄} = ht_M̄ + ht_N̄`

★★引き戻し `x_F^* M̄` は `§9-742`〜`§9-745`（計量の引き戻しと `APic` の群準同型）で
無条件に入っている。★★★`deg_F` も `§9-782` で無条件になった。
したがって**本ファイルの主張はすべて無条件である**。

## ★残っている段（明示）

★★`X(ℚ̄)` の上の関数として述べるには、底変換

    `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)`（計量版）

が要る。★★★因子の側ではそれは `degNormalized_baseChange` として証明済みであり、
計量版は「有限部分もアルキメデス部分もちょうど `[K:F]` 倍になる」ことに帰着する
——正規化 `/[F:ℚ]` がそれを打ち消す。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField
open ABC3.Found.GenEll

/-- ★★★★★★★★★★**`ht_M̄(x_F) = deg_F(x_F^* M̄)`（無条件）**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`MetricHeight.lean` の `htMetric` と違い、[Szp] Prop 1.1 を仮定しない。 -/
noncomputable def htMetricU {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (M : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) : ℝ :=
  degAPicM F (APicMPullback xF (APicM.mk M))

/-- ★★★★★★★★★★**高さは加法的である**（無条件）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★機構は `§9-745` の `APicMPullback` が群準同型であることと、
`§9-782` の `degAPicM_mul` の合成だけである。 -/
theorem htMetricU_mul {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (M N : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F (M.mul N) xF = htMetricU F M xF + htMetricU F N xF := by
  show degAPicM F (APicMPullback xF (APicM.mk (M.mul N))) = _
  have h : APicM.mk (M.mul N) = APicM.mk M * APicM.mk N := rfl
  rw [h, map_mul, degAPicM_mul]
  rfl

/-- ★★**自明な算術直線束の高さは `0`**（非空虚性の witness、無条件）。 -/
@[simp] theorem htMetricU_one {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F (AInv.one X) xF = 0 := by
  show degAPicM F (APicMPullback xF (APicM.mk (AInv.one X))) = _
  have h : APicM.mk (AInv.one X) = (1 : APicM X) := rfl
  rw [h, map_one, degAPicM_one]

/-- ★★★**高さは等長同型類だけで決まる**（無条件）。 -/
theorem htMetricU_congr {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    {M N : AInv X} (h : Isometric M.carrier N.carrier)
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F M xF = htMetricU F N xF := by
  show degAPicM F (APicMPullback xF (APicM.mk M)) = _
  rw [APicM.mk_eq_mk h]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def htMetricU.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(ht_M̄(x) = deg_F(x_F^* M̄)——[Szp] を仮定しない形)",
    sectionId := "genell-def-1-1-ii" }

def htMetricU_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(高さの加法性——[Szp] を仮定しない形)",
    sectionId := "genell-def-1-1-ii" }

def htMetricU_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "APicMPullback(引き戻しが APic の群準同型、§9-745)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicMPullback") 4,
    .citation "[ABC3]" "degAPicM_mul(deg_F は準同型、§9-782)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degAPicM_mul") 4,
    .implicitStep
      ("★★X(ℚ̄) の上の関数として述べるには底変換 deg_K(L|_{Spec 𝓞_K}) = deg_F(L) の" ++
       "計量版が要る。★★★因子の側では degNormalized_baseChange として証明済みであり、" ++
       "計量版は「有限部分もアルキメデス部分もちょうど [K:F] 倍になる」ことに帰着する" ++
       "——正規化 /[F:ℚ] がそれを打ち消す") 4 ]

end ABC3.Found.Arakelov
