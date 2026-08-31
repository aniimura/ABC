/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackPic
import ABC3.Found.GenEll.ProductFormula
import ABC3.Meta.Claim

/-!
# **計量表示の高さ** `ht_M̄(x) = deg_F(x_F^* M̄)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★原文の構造をそのまま写す

`Definition 1.1 (ii)` の `deg_F` は**因子の側で定義される**:

    `deg_F : ADiv(F)/APrc(F) → ℝ`,  `v ↦ log(q_v)`（非アルキメデス）, `v ↦ 1`（アルキメデス）

そして原文は

> as is well-known, there is a natural isomorphism
> `ADiv(F)/APrc(F) ≅ APic(Spec(O_F))` — cf. [Szp], Proposition 1.1

と**引用**して、それを `APic(Spec 𝒪_F) → ℝ` と読み替える。

★★**原文が証明せず引用している段はここだけである**。したがって本ファイルは
その同型を**仮定（`SzpData`）として担ぐ**——それが原文の依存関係の忠実な写しである。

## ★在庫との関係

| もの | 場所 | 状態 |
|---|---|---|
| `ADiv` / `APrc` / `deg` / 正規化 | `ArithDiv.lean` | ✅ 無条件 |
| `deg_F : ADiv/APrc →+ ℝ`（積公式で well-defined） | `ProductFormula.lean` | ✅ 無条件 |
| `deg_K(L|_K) = deg_F(L)` | `BaseChange.lean` | ✅ 無条件 |
| `φ^* : APicM Y →* APicM X` | `PullbackPic.lean`（§9-745） | ✅ 無条件 |
| **`ADiv/APrc ≅ APic(Spec 𝒪_F)`** | ★**[Szp] Prop 1.1** | ★仮定 |

## ★★★★★★★★出るもの

    `ht_M̄(x_F) ≝ deg_F(x_F^* M̄)`,  `ht_{M̄ ⊗ N̄} = ht_M̄ + ht_N̄`

★加法性は **`§9-745` の `APicMPullback` が群準同型であること**と
`degNormalizedAPic` が加法準同型であることの合成であり、**それ以上は要らない**。
★★これが `§9-742`〜`§9-745` を積み上げた目的である。

## ★残っている段（明示）

★★`X(ℚ̄)` の上に降ろすには、`SzpData` の**族**が底変換と両立すること
（原文の `deg_K(L|_{Spec 𝒪_K}) = deg_F(L)` の計量版）が要る。
★因子の側ではそれは `degNormalized_baseChange` として**証明済み**なので、
残るのは「[Szp] の同型の族が底変換と両立する」ことだけである。
★★★`AMetricPullback` の**関手性**（`(f ≫ g)^* ≅ f^* ∘ g^*`）も、その段で要る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★[Szp] Prop 1.1 —— 原文自身の引用 -/

/-- ★★★★★★**[Szp] Prop 1.1 の同型**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★原文は `ADiv(F)/APrc(F) ≅ APic(Spec(O_F))` を
「as is well-known … cf. [Szp], Proposition 1.1」と**引用**する。
★★したがってこれは依存グラフの**辺**であって、我々が証明すべき節点ではない。
本構造体はその引用を型で担ぐ。

★★★注意: `APicM` は**前層の水準**である（`AMetricPic.lean` の逸脱の記録）。 -/
structure SzpData (F : Type) [Field F] [NumberField F] where
  /-- `APic(Spec 𝒪_F) ≅ ADiv(F)/APrc(F)`。 -/
  iso : APicM (Spec (CommRingCat.of (𝓞 F))) ≃* Multiplicative (APicOF F)

/-! ## ★★★★★★★`deg_F` を `APic` の上で読む -/

/-- ★★★★★★★**`deg_F : APic(Spec 𝒪_F) → ℝ`**（正規化した版）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★因子の側の `degNormalizedAPic` を [Szp] の同型で運んだものである。 -/
noncomputable def degMetric (F : Type) [Field F] [NumberField F] (szp : SzpData F)
    (L : APicM (Spec (CommRingCat.of (𝓞 F)))) : ℝ :=
  degNormalizedAPic (Multiplicative.toAdd (szp.iso L))

/-- ★★★★**`deg_F` は準同型である**。 -/
theorem degMetric_mul (F : Type) [Field F] [NumberField F] (szp : SzpData F)
    (L M : APicM (Spec (CommRingCat.of (𝓞 F)))) :
    degMetric F szp (L * M) = degMetric F szp L + degMetric F szp M := by
  show degNormalizedAPic (Multiplicative.toAdd (szp.iso (L * M))) = _
  rw [map_mul]
  exact map_add (degNormalizedAPic (F := F)) _ _

@[simp] theorem degMetric_one (F : Type) [Field F] [NumberField F] (szp : SzpData F) :
    degMetric F szp 1 = 0 := by
  show degNormalizedAPic (Multiplicative.toAdd (szp.iso 1)) = _
  rw [map_one]
  exact map_zero (degNormalizedAPic (F := F))

/-! ## ★★★★★★★★★★計量表示の高さ -/

/-- ★★★★★★★★★★**`ht_M̄(x) = deg_F(x_F^* M̄)`**（計量表示）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★引き戻し `x_F^* M̄` は `§9-742`（計量の引き戻し）と `§9-745`
（`APic` の群準同型）で無条件に入っている。
★★[Szp] の同型だけが仮定である。 -/
noncomputable def htMetric {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (szp : SzpData F) (M : AInv X)
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) : ℝ :=
  degMetric F szp (APicMPullback xF (APicM.mk M))

/-- ★★★★★★★★★★**高さは加法的である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★機構は `§9-745` の `APicMPullback` が**群準同型**であることと、
`degNormalizedAPic` が加法準同型であることの合成だけである。
★★★これが `§9-742`〜`§9-745` を積み上げた目的である
——引き戻しの乗法性（`§9-743`）と単位の等長性（`§9-745`）が
ここで初めて**高さの加法性**として現れる。 -/
theorem htMetric_mul {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (szp : SzpData F) (M N : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetric F szp (M.mul N) xF = htMetric F szp M xF + htMetric F szp N xF := by
  show degMetric F szp (APicMPullback xF (APicM.mk (M.mul N))) = _
  have h : APicM.mk (M.mul N) = APicM.mk M * APicM.mk N := rfl
  rw [h, map_mul, degMetric_mul]
  rfl

/-- ★★**自明な算術直線束の高さは `0`**（非空虚性の witness）。 -/
@[simp] theorem htMetric_one {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (szp : SzpData F) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetric F szp (AInv.one X) xF = 0 := by
  show degMetric F szp (APicMPullback xF (APicM.mk (AInv.one X))) = _
  have h : APicM.mk (AInv.one X) = (1 : APicM X) := rfl
  rw [h, map_one, degMetric_one]

/-- ★★★**高さは等長同型類だけで決まる**。

★原文の `ht` が `APic` の上の関数であること（束そのものではなく類の関数）に当たる。 -/
theorem htMetric_congr {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (szp : SzpData F) {M N : AInv X} (h : Isometric M.carrier N.carrier)
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetric F szp M xF = htMetric F szp N xF := by
  show degMetric F szp (APicMPullback xF (APicM.mk M)) = _
  rw [APicM.mk_eq_mk h]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def SzpData.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)([Szp] Prop 1.1 の同型を仮定として担ぐ形)",
    sectionId := "genell-def-1-1-ii" }

def degMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F を APic(Spec 𝓞_F) の上で読む——[Szp] を仮定)",
    sectionId := "genell-def-1-1-ii" }

def htMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(ht_M̄(x) = deg_F(x_F^* M̄)——[Szp] を仮定)",
    sectionId := "genell-def-1-1-ii" }

def htMetric_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(高さの加法性——[Szp] を仮定)",
    sectionId := "genell-def-1-1-ii" }

def htMetric_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[Szp]" "Proposition 1.1(ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F))"
      (.absent ("2026-08-28 実測: mathlib(lean/.lake/packages/mathlib/Mathlib/)を " ++
        "arakelov(大小無視)で 0 ファイル、" ++
        "ArithmeticDivisor / arithmeticPicard / ArakelovDivisor で 0 件。" ++
        "★算術 Picard 群(計量つき)そのものが mathlib に無い")) 4,
    .citation "[ABC3]" "APicMPullback(引き戻しが APic の群準同型、§9-745)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicMPullback") 4,
    .citation "[ABC3]" "degNormalizedAPic(deg_F が ADiv/APrc の上で well-defined——積公式)"
      (.inProject "ABC3" "ABC3.Found.GenEll.degNormalizedAPic") 4,
    .implicitStep
      ("★原文は ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F) を " ++
       "「as is well-known … cf. [Szp], Proposition 1.1」と**引用**する。" ++
       "したがってこれは依存グラフの**辺**であって、我々が証明すべき節点ではない。" ++
       "★★本ファイルはその引用を型(SzpData)で担ぐ") 4,
    .implicitStep
      ("★★★残っている段: X(ℚ̄) の上に降ろすには SzpData の**族**が底変換と" ++
       "両立すること(原文の deg_K(L|_K) = deg_F(L) の計量版)が要る。" ++
       "因子の側ではそれは degNormalized_baseChange として証明済み。" ++
       "★AMetricPullback の関手性((f ≫ g)^* ≅ f^* ∘ g^*)もその段で要る") 4 ]

end ABC3.Found.Arakelov
