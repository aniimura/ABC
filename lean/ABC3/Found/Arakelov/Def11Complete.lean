/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.HeightBaseChange
import ABC3.Found.Arakelov.ContMetric
import ABC3.Found.GenEll.AlgPointClass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★[GenEll] Definition 1.1 —— **項目全体が揃った**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3–4。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★本ファイルが足す最後の 1 本 —— `ht_M̄ : X(ℚ̄) → ℝ`

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

`§9-799` で `ht_M̄` が `x_F` の取り方に依らないことを示した（`htMetricU_baseChange`）。
★本ファイルはそれを使って `X(ℚ̄) = ⋃_F X(F)`（`AlgPointAnyClass`）の上へ**降ろす**
——`Quot.lift` 一発である。★★これで原文の

    `ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ`

が **`X(ℚ̄)` 上の関数として** well-defined になる。

## ★★★★★★★★★★主張 —— 原文の条と宣言の対応

### (i)

| 原文 | 宣言 |
|---|---|
| `X^arc`（複素点・位相・コンパクト性） | `complexPoints` / `arcTopology` / `compactSpace_arc` |
| `ι_X`（複素共役、対合） | `conjPoint` / `conjPoint_conjPoint` |
| 直線束 `L` on `X` | `AMetric.sheaf` ＋ `AMetric.triv`（`IsLocallyTrivial`） |
| `L^arc` 上の hermitian 計量 `|−|_L` | `AMetric.metric`（`LocalMetric`）／`AMetric.norm` |
| ——`ι_X` と**両立** | ★`AMetric.IsConjCompatible` ／ `AMetric.norm_conjPoint`（`§9-800`） |
| ——**連続**（`X^arc` は解析空間） | ★`AMetric.IsContinuous` ／ `AMetric.continuous_norm`（`§9-802`） |
| 射 `L̄ → M̄` | `IsAHom` ／ `AHom` |
| `Γ(L̄)` は `Ō_X → L̄` の射の集合 | `AMetric.gamma` ／ `isAHom_one_iff` |
| テンソル積 | `AMetric.mul`（`LocalMetric.tensor`＋`IsTensorOf`）／`AMetric.normOf_mul` |
| 同型類が群 `APic(X)` をなす | ★`APicA X = APicC X ⊓ APicK X`（`§9-801`/`§9-802`） |

### (ii)

| 原文 | 宣言 |
|---|---|
| 引き戻し `φ^*L̄ = L̄|_Y` | `AMetricPullback` ／ `AInv.pullback` ／ `APicMPullback` |
| ——共役両立・連続を保つ | `AMetricPullback_isConjCompatible` ／ `AMetricPullback_isContinuous` |
| `V(F) = V(F)^non ⊔ V(F)^arc` | `FinitePlace` ／ `InfinitePlace` |
| `q_v`・`ord_v`・`[F_v:ℝ]` | `residueCard` ／ `ordv` ／ `InfinitePlace.mult` |
| 算術因子 `ADiv(F)`・effective | `ADiv` ／ `ADiv.IsEffective` |
| `ADiv(f)`・`APrc(F)` | `principalADiv` ／ `APrc` |
| `ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)` | ★`SzpData`（**[Szp] Prop 1.1 の引用を型で受ける**） |
| `deg_F`（`v ↦ log q_v`／`v ↦ 1`） | `deg` ／ `degHom` |
| 正規化 `deg_F ≝ (1/[F:ℚ])·deg_F` | `degNormalized` |
| `deg_F : APic(Spec 𝓞_F) → ℝ`（★**[Szp] を経由しない**） | ★`degAPicM` ／ `degAPicM_mul`（`§9-782`） |
| `deg_K(L̄|_{Spec 𝓞_K}) = deg_F(L̄)` | ★`degAPicM_baseChange`（`§9-799`） |
| `X(ℚ̄) = ⋃_{[F:ℚ]<∞} X(F)` | `AlgPointAnyClass` |
| `ht_M̄(x) ≝ deg_F(x_F^* M̄)`、`x_F` に依らない | ★`htMetricAlg`（**本ファイル**）／`htMetricU_baseChange` |

## ★★★★★逸脱の記録

### 1. 前層の水準で語る

`AMetric` の台は `X.PresheafOfModules` ＋ `IsLocallyTrivial` であって、
連接層の水準ではない（`AMetricPic.lean` の記録）。★`Spec 𝓞_F` 上では
`§9-779` の測定（`Γ_pre(L)` は可逆 `Γ(X,⊤)`-加群）でそのまま次数が取れる。

### 2. 射の条件を**強い側**で取った

原文は「**局所的に** `|−|_L ≤ 1` の切断が `|−|_M ≤ 1` へ移る」だが、
`IsAHom` は `|φ(s)|_M ≤ |s|_L`（作用素ノルム `≤ 1`）を要求する。
★直線束の上では同値だが**その同値は証明していない**。向きは強い側なので下流は弱くならない。

### 3. `X` に正規性・`ℤ`-固有性・`ℤ`-平坦性を課していない

それらは §1 の地の文の標準仮定であって、定義そのものには要らない。
★要る所（`X^arc` のコンパクト性）では `compactSpace_arc` が固有性から出す。

### 4. `X^arc` は「複素点＋`arcTopology`」である

原文の `X^arc` は**コンパクト正規複素解析空間**だが、
★mathlib に解析化（analytification / GAGA / complex analytic space）は**無い**
（2026-08-16 実測、いずれも 0 件）。★★本実装の `X^arc` は
`complexPoints X`（＝`Spec ℂ ⟶ X`）に `arcTopology` を入れたものであり、
計量の連続性はこの位相についてである。

### 5. [Szp] Prop 1.1 は**引用として受け、使っていない**

原文自身が「as is well-known … cf. [Szp], Proposition 1.1」と引用する。
★`SzpData` がそれを型で担ぐ。★★しかし本実装の `deg_F`（`degAPicM`）は
**その同型を経由せずに**直接作ってある（`§9-776`〜`§9-782`）ので、
`ht_M̄` は `SzpData` を**仮定しない**。

### 6. `deg_F`・`ht_M̄` は `APicM` 全体の上で定義してある

原文の `APic(X)` は共役両立かつ連続な計量を持つ類（`APicA X`）だが、
★`degAPicM` / `htMetricAlg` はその制限を課さずに定義できる（**より広い**）。
★★`APicA X` へ制限しても同じ値である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★★★★★`ht_M̄` を `X(ℚ̄)` の上へ降ろす -/

/-- ★★定義体つきの点における高さ。 -/
noncomputable def htMetricAny {X : Scheme.{0}} (M : AInv X) (p : AlgPointAny X) : ℝ :=
  @htMetricU X p.fld p.instField p.instNF M p.map

/-- ★★★★★★★★**定義体を上げても値は変わらない**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★これが `§9-799` の `htMetricU_baseChange` そのものである。 -/
theorem htMetricAny_respects {X : Scheme.{0}} (M : AInv X) {p q : AlgPointAny X}
    (h : AlgPointAnyRel p q) : htMetricAny M p = htMetricAny M q := by
  cases h with
  | base F K xF => exact (htMetricU_baseChange F K M xF).symm

/-- ★★★★★★★★★★**`ht_M̄ : X(ℚ̄) → ℝ`**（無条件）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ`

★★これが原文の高さである。`x_F` の取り方に依らないので `X(ℚ̄)` の上の関数になる。
★★★[Szp] Prop 1.1 を**仮定していない**。 -/
noncomputable def htMetricAlg {X : Scheme.{0}} (M : AInv X) : AlgPointAnyClass X → ℝ :=
  Quot.lift (htMetricAny M) (fun _ _ h => htMetricAny_respects M h)

@[simp] theorem htMetricAlg_mk {X : Scheme.{0}} (M : AInv X) (F : Type) [Field F]
    [NumberField F] (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricAlg M (AlgPointAnyClass.mk (algPointAny F xF)) = htMetricU F M xF := rfl

/-- ★★★★★★**高さは加法的である**（`X(ℚ̄)` の上で）。 -/
theorem htMetricAlg_mul {X : Scheme.{0}} (M N : AInv X) (x : AlgPointAnyClass X) :
    htMetricAlg (M.mul N) x = htMetricAlg M x + htMetricAlg N x := by
  induction x using AlgPointAnyClass.ind with
  | _ p => exact @htMetricU_mul X p.fld p.instField p.instNF M N p.map

/-- ★★自明な算術直線束の高さは `0`（非空虚性の witness）。 -/
@[simp] theorem htMetricAlg_one {X : Scheme.{0}} (x : AlgPointAnyClass X) :
    htMetricAlg (AInv.one X) x = 0 := by
  induction x using AlgPointAnyClass.ind with
  | _ p => exact @htMetricU_one X p.fld p.instField p.instNF p.map

/-- ★★★★**高さは等長同型類だけで決まる**。 -/
theorem htMetricAlg_congr {X : Scheme.{0}} {M N : AInv X}
    (h : Isometric M.carrier N.carrier) (x : AlgPointAnyClass X) :
    htMetricAlg M x = htMetricAlg N x := by
  induction x using AlgPointAnyClass.ind with
  | _ p => exact @htMetricU_congr X p.fld p.instField p.instNF M N h p.map

/-! ## ★出典の紐付け(`.src`) -/

def htMetricAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(ht_M̄ : X(ℚ̄) → ℝ——x_F の取り方に依らない形、[Szp] 不使用)",
    sectionId := "genell-def-1-1-ii" }

/-- ★★★★★★★★★★**[GenEll] Definition 1.1**（算術直線束・`APic(X)`・`deg_F`・`ht`）
—— 実装された。

原文 (GenEll p.3):
> The isomorphism classes of arithmetic line bundles on X, together with the operation
> of tensor product, thus determine a group APic(X).

★対応表と逸脱の記録は本ファイル冒頭の docstring にある。

★★**過去の取り下げとの違い**（2026-08-28）: 以前 `ArithPic` の設計で
本項目全体の `.src` を置き、取り下げた（`Definition11.lean` の記録）。
理由は「群法則が**テンソル積ではなかった**」ことである
——`TorsorMetric.tensor` は Green 関数を足すだけで、基準計量が
`Classical.choice` だったため 2-コサイクルの差が残っていた。

★★★**いまの設計ではその差は無い**: `AMetric.mul` は `LocalMetric.tensor` であり、
`isTensorOf_tensor` が `h_{A⊗B} = h_A · h_B` を、`AMetric.normOf_mul` が
`|s ⊗ t| = |s|·|t|` を**定理として**与える。群法則は本物のテンソル積である。 -/
def definition_1_1_metric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1",
    sectionId := "genell-def-1-1-i" }

def definition_1_1_metric.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[Szp]" "Proposition 1.1(ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)——原文自身が引用)"
      (.absent ("[Szp] は ResearchPaper/papers.json に登記済みだが未目視。" ++
        "本実装は SzpData で型として受けるだけで、deg_F / ht には使っていない")) 4,
    .citation "[ABC3]" "APicA(共役両立かつ連続な類がなす群 APic(X)、§9-801/§9-802)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicA") 3,
    .citation "[ABC3]" "AMetric.norm_conjPoint(|s|(ι_X p) = |s|(p)、§9-800)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_conjPoint") 3,
    .citation "[ABC3]" "AMetric.continuous_norm(|s|_L̄ は X^arc 上連続、§9-802)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.continuous_norm") 3,
    .citation "[ABC3]" "AMetric.normOf_mul(|s ⊗ t| = |s|·|t|——群法則が本物のテンソル積)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.normOf_mul") 3,
    .citation "[ABC3]" "degAPicM / degAPicM_mul(deg_F : APic(Spec 𝓞_F) → ℝ、§9-782)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degAPicM") 4,
    .citation "[ABC3]" "degAPicM_baseChange(deg_K(L|Spec 𝓞_K) = deg_F(L)、§9-799)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degAPicM_baseChange") 4,
    .citation "[ABC3]" "ADiv / APrc / principalADiv / deg / degNormalized(算術因子の側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ADiv") 4,
    .citation "[ABC3]" "AlgPointAnyClass(X(ℚ̄) = ⋃_F X(F))"
      (.inProject "ABC3" "ABC3.Found.GenEll.AlgPointAnyClass") 4,
    .implicitStep
      ("★逸脱 4: 原文の X^arc は**コンパクト正規複素解析空間**だが、mathlib に " ++
       "analytification / GAGA / complex analytic space は無い(2026-08-16 実測、0 件)。" ++
       "★★本実装の X^arc は complexPoints X に arcTopology を入れたものであり、" ++
       "計量の連続性はこの位相についてである") 3,
    .implicitStep
      ("★逸脱 2: 射の条件を強い側(作用素ノルム ≤ 1)で取った。" ++
       "原文の「局所的に |−|_L ≤ 1 の切断が |−|_M ≤ 1 へ移る」との同値は" ++
       "直線束の上では成り立つが**証明していない**。向きは強い側なので下流は弱くならない") 3 ]

end ABC3.Found.Arakelov
