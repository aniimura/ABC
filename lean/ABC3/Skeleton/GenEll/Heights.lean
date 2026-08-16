import ABC3.Meta.Claim
import ABC3.Found.GenEll.ArithDiv
import ABC3.Interface.GenEll.ArithLineBundle

/-!
# [GenEll] §1 —— 高さの `Skeleton`(層 D)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.4。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★層の現在地

| 層 | 中身 | 場所 | 状態 |
|---|---|---|---|
| A | `ADiv(F)` / `deg_F` / `ord_v` / 主因子 `ADiv(f)` / `APrc(F)` | `Found/GenEll/ArithDiv.lean` | ★**実装済**(`sorry` 無し) |
| B・C | 算術直線束と引き戻し(スキーム上の直線束・解析化・hermitian 計量) | `Interface/GenEll/ArithLineBundle.lean` | `waiting` |
| **D** | **`ht` の定義と、`F` の取り方に依らないこと** | ★**本ファイル** | 下記 |

## ★★★`sorry` を消した方法と、その正直な意味(必読)

本ファイルの `sorry` は **0** だが、それは「証明した」という意味ではない。

★**内容は `Interface` の仮説として輸入した。** 具体的には
`PulledBackClassData` が `base_change_invariant`(不変性)と `height_eq`(高さの一致)を
**posit** しており、本ファイルの定理はそれを取り出しているだけである。

★これは `tools/check.mjs` 冒頭 **B5 が名指しする穴**そのものである——
> 「`Interface` に条件を posit すれば sorry は消える(実例: IUTchIII cor_3_12、2026-08-14)」

`Skeleton/IUTchIII/Cor312.lean` の `cor_3_12_inequality` が同じ形であり、
そこでは `.needs` に `.implicitStep "… Interface の仮説として輸入した(我々は証明していない)"`
と明記している。**本ファイルも同じ規律に従う**——各定理の `.needs` に明記した。

★**したがって「`sorry` が 0 になった」を進捗と読んではならない。**
本当の進捗は `Interface` が `waiting` でなくなったときに起きる。
`node tools/check.mjs` の「Interface 実装待ち」の行がそれを見せる。

## ★仮説を強めていないことの確認(`PLAN.md` A6)

`Interface` の 2 つの仮説は、どちらも
`DefinedOver F x`(`x` が `F` 上で定義される)と `[Algebra F K]` の下でだけ主張されている。
★無条件の等式として posit すると原文より強くなるので、そうしていない。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll ABC3.Interface.GenEll NumberField

section Height

variable (D : PulledBackClassData)

/-- **高さ** `ht_M̄`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

原文の定義は `ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ`、`x ∈ X(ℚ̄) = ∪_{[F:ℚ]<∞} X(F)`。

★**`Interface` が `height` として posit したものを取り出しているだけ**である。
中身(層 B・C)を作ったわけではない。 -/
noncomputable def ht (M : D.Bundle) (x : D.Point) : ℝ := D.height M x

/-- ★★**高さが `F` の取り方に依らない**こと。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

原文は `deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`(正規化した `deg`)と述べ、
それによって `ht` が `X(ℚ̄)` の上で定まる。

★**証明していない**——`Interface` の `base_change_invariant` を取り出しているだけである
(上の docstring と `.needs` を参照)。 -/
theorem degNormalized_base_change
    (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
    (M : D.Bundle) (x : D.Point) (hF : D.DefinedOver F x) (hK : D.DefinedOver K x) :
    D.degOver K M x = D.degOver F M x :=
  D.base_change_invariant F K M x hF hK

/-- ★`ht` が、`x` を定義するどの数体の上での次数とも一致すること。

★これが「`ht` が well-defined である」ことの内容である。
**証明していない**——`Interface` の `height_eq` を取り出しているだけ。 -/
theorem ht_eq_degOver (F : Type) [Field F] [NumberField F]
    (M : D.Bundle) (x : D.Point) (hF : D.DefinedOver F x) :
    ht D M x = D.degOver F M x :=
  D.height_eq F M x hF

/-- ★**2 つの数体の上での次数が一致する**——`ht` を経由した形。

`ht_eq_degOver` を 2 回使えば出る。★`Interface` の 2 仮説が
互いに整合していることの確認でもある(片方だけでは出ない)。 -/
theorem degOver_eq_of_definedOver
    (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    (M : D.Bundle) (x : D.Point) (hF : D.DefinedOver F x) (hK : D.DefinedOver K x) :
    D.degOver F M x = D.degOver K M x := by
  rw [← ht_eq_degOver D F M x hF, ← ht_eq_degOver D K M x hK]

end Height

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def degNormalized_base_change.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

/-- ★★**この定理は証明していない**——内容を `Interface` の仮説として輸入した。 -/
def degNormalized_base_change.needs : List ProofObligation :=
  [ .implicitStep
      "★**内容を Interface の仮説として輸入した(我々は証明していない)**。`PulledBackClassData.base_change_invariant` がそれである。check.mjs 冒頭 B5 が名指しする穴と同じ形で、`Skeleton/IUTchIII/Cor312.lean` の `cor_3_12_inequality` と同型" 4,
    .citation "[Szp]" "§1.1(算術因子と次数の基本性質)"
      (.absent "mathlib 全体を `Arakelov` / `arithmetic.*line bundle` / `LineBundle` / `analytification` で grep、いずれも 0 件(2026-08-16)") 4,
    .derivation
      "本当に要るのは素点の分岐・剰余次数の関係 `Σ_{w|v} e_w f_w = [K:F]`。mathlib に `Ideal.sum_ramification_inertia` 系があるはずだが**未測定**" 4 ]

def ht.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(高さ ht_M̄)",
    sectionId := "genell-ht" }

/-- ★★**`ht` は構成していない**——`Interface` が posit した `height` を取り出しただけ。 -/
def ht.needs : List ProofObligation :=
  [ .implicitStep
      "★**`ht` の構成そのものを Interface に posit させた**。本当に構成するには、`x_F^* M̄ ∈ APic(Spec 𝒪_F)` を作る操作(層 B・C)が要る" 4,
    .otherPaper "[GenEll]" "Definition 1.1, (i)(算術直線束——層 B・C)" 3,
    .otherPaper "[GenEll]" "Definition 1.1, (ii)(引き戻し)" 3 ]

def ht_eq_degOver.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(高さ ht_M̄)",
    sectionId := "genell-ht" }

def ht_eq_degOver.needs : List ProofObligation :=
  [ .implicitStep
      "★**内容を Interface の仮説として輸入した(`height_eq`)**。原文は `ht` の定義の直後に「any morphism that gives rise to x」と書くだけで、well-defined 性を証明していない" 4 ]

def degOver_eq_of_definedOver.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

/-- ★これは**本当に証明した**唯一の定理である(`Interface` の 2 仮説から導いた)。 -/
def degOver_eq_of_definedOver.needs : List ProofObligation :=
  [ .derivation "`ht_eq_degOver` を 2 回。★Interface の仮説からの導出であり、ここは我々が書いた" 4 ]

end ABC3.Skeleton.GenEll
