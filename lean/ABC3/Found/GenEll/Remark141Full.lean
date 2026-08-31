/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Remark141
import ABC3.Found.GenEll.ZModelOfProjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★[GenEll] Remark 1.4.1 —— 3 条を 1 本に（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★★★★★★★★★★★★★★★内訳

原文の Remark は 2 文である:

1. 高さ関数の **BD-class の理論は `X_ℚ` だけに依る**
2. したがって理論は **`ℚ` 上の正規射影的スキーム `Y`** に適用できる

★★本ファイルはその両方を 1 本にまとめる:

| 条 | 内容 | 実装 |
|---|---|---|
| (1) | 同じモデルの上——類が同じなら BD-同値 | `htArith_bdeq_of_pullbackAPic_eq`（`HeightClass.lean`） |
| (2) | **2 つのモデルにまたがっても**——点の対応を挟めば同じ | `htArith_bdeq_of_pullbackAPic_eq'`（`Remark141.lean`） |
| (3) | 射影的な `Y` には **`Spec ℤ` 上固有な `ℤ`-モデル**がある | ★`§9-953`（`ZModelOfProjective.lean`、本日） |

★★★★★★**(3) が本日（2026-08-29）取れたので、第 2 文の側が揃った**。

## ★★★★★逸脱の記録（`CLAUDE.md` の「逸脱」）

### 1. 「だけに依る」の中身

原文は「`X_ℚ` だけに依る」。★本実装は **「引き戻した算術因子の類だけに依る」**である。
★★2 つの `ℤ`-モデルのずれ（＝**垂直因子**の差）は `Remark141.lean` の梯子で扱っている:

| 版 | 仮定 | 定数 |
|---|---|---|
| `remark_1_4_1_of_descent` | 類が一致 | `0` |
| `remark_1_4_1_of_descent_twist` | 差が `Spec ℤ` から来る | `deg_ℚ(N)` |
| `remark_1_4_1_of_descent_bound` | 差のイデアルが `n` を含む | `log n` |
| `remark_1_4_1_of_descent_comparable` | `n` で互いに比較できる | `log n` |
| `htArith_bdeq_of_chartwise_localization` | ★**一様な `n` をスキームから取る**（第 365） | `log n` |

### 2. 定数

原文は BD-class（定数差を許す）。★構成の側では **`C = 0`**——原文より強い。
機構は `deg` が主算術因子の上で消えること（積公式）。

### 3. (3) で受けているもの（★明示）

★★★**`Y ⟶ ℙᴺ_ℤ` は与えられたものとしている**。
原文の「`ℚ` 上正規射影的」から `Y ↪ ℙᴺ_ℚ` を得て、さらに `ℙᴺ_ℚ ⟶ ℙᴺ_ℤ`（底変換の射影）を
合成する段は**含めていない**——後者は `Proj` の底変換の射を作る段であり、
本プロジェクトにまだ無い。
★**そこから先（像が `ℤ`-モデルで、固有・分離・コンパクトであること）は本ファイルが持つ**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★★★★★★★★★★★★★★★★★★★★3 条 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Remark 1.4.1** —— 3 条。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

* **(1)** 同じモデルの上——引き戻した類が同じなら高さは BD-同値（定数は `0`）
* **(2)** **2 つのモデルにまたがっても**——点の対応 `ePt` を挟めば同じ
* **(3)** 射影的な `Y` には **`Spec ℤ` 上固有な `ℤ`-モデル**がある（第 2 文の側）

★★★★★逸脱はファイル冒頭に記録した。 -/
theorem remark_1_4_1_full :
    -- ★(1) 同じモデルの上: 引き戻した類が同じなら高さは BD-同値
    (∀ {X : Scheme.{0}} (D E : ArithCartier X) (F : Type) [Field F] [NumberField F],
        (∀ xF : specRingOfIntegers F ⟶ X, pullbackAPic F D xF = pullbackAPic F E xF) →
        BDeq (fun xF => htArith F D xF) (fun xF => htArith F E xF))
    -- ★★(2) 2 つのモデルにまたがっても、点の対応を挟めば同じ
  ∧ (∀ {X X' : Scheme.{0}} (D : ArithCartier X) (E : ArithCartier X')
        (F : Type) [Field F] [NumberField F]
        (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X')),
        (∀ xF, pullbackAPic F D xF = pullbackAPic F E (ePt xF)) →
        BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF)))
    -- ★★★(3) 射影的な Y には Spec ℤ 上固有な ℤ-モデルがある
  ∧ (∀ {Y : Scheme.{0}} {N : ℕ}
        (e : Y ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)),
        toZModel e ≫ zModelι e = e ∧
        IsProper (zModelι e ≫
          Proj.toSpecZero (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) ∧
        Scheme.IsSeparated (zModel e) ∧ CompactSpace (zModel e : Type)) := by
  refine ⟨?_, ?_, ?_⟩
  · intro X D E F _ _ h
    exact htArith_bdeq_of_pullbackAPic_eq F D E h
  · intro X X' D E F _ _ ePt h
    exact htArith_bdeq_of_pullbackAPic_eq' F D E ePt h
  · intro Y N e
    exact ⟨toZModel_zModelι e, isProper_zModel e, isSeparated_zModel e, compactSpace_zModel e⟩

/-! ## ★出典の紐付け(`.src`) -/

def remark_1_4_1_full.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Remark 1.4.1",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1_full.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_bdeq_of_pullbackAPic_eq((1) 同じモデルの上)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_pullbackAPic_eq") 3,
    .citation "[ABC3]" "htArith_bdeq_of_pullbackAPic_eq'((2) 2 つのモデルにまたがって)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_pullbackAPic_eq'") 3,
    .citation "[ABC3]" "isProper_zModel((3) 射影埋め込みから ℤ-モデル、§9-953)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isProper_zModel") 4,
    .citation "[ABC3]"
      "htArith_bdeq_of_chartwise_localization(一様な n をスキームから取る段、第 365)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_chartwise_localization") 4,
    .implicitStep
      ("★★★★★逸脱 1(「だけに依る」の中身): 原文は「X_ℚ だけに依る」だが、" ++
       "本実装は「引き戻した算術因子の**類**だけに依る」である。" ++
       "2 つの ℤ-モデルのずれ(＝垂直因子の差)は Remark141.lean の梯子" ++
       "(descent / twist / bound / comparable)と ChartwiseUniform.lean" ++
       "(一様な n をスキームから取る段、第 365)で扱っている") 6,
    .implicitStep
      ("★★逸脱 2(定数): 原文は BD-class(定数差を許す)だが、構成の側では C = 0" ++
       "——原文より強い。機構は deg が主算術因子の上で消えること(積公式)である") 3,
    .implicitStep
      ("★★★★逸脱 3((3) で受けているもの): **Y ⟶ ℙᴺ_ℤ は与えられたものとしている**。" ++
       "原文の「ℚ 上正規射影的」から Y ↪ ℙᴺ_ℚ を得て ℙᴺ_ℚ ⟶ ℙᴺ_ℤ(底変換の射影)を" ++
       "合成する段は含めていない——Proj の底変換の射を作る段は本プロジェクトにまだ無い。" ++
       "★そこから先(像が ℤ-モデルで、固有・分離・コンパクトであること)は本ファイルが持つ") 5 ]

end ABC3.Found.GenEll
