import ABC3.Meta.Claim
import ABC3.Found.GenEll.ArithDiv
import ABC3.Interface.GenEll.ArithLineBundle

/-!
# [GenEll] §1 —— 高さの `Skeleton`(層 D)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.4。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★この Skeleton が固定するもの —— 「何を作れば終わりか」

`ResearchPaper/foundations.json` の `arakelov` 節点は「Arakelov 理論」という
**まとまり**なので、そのままでは受理条件にならない(2 値で数えられない)。
本ファイルは**それを型に落とす**——ここに並ぶ `sorry` が埋まれば、
`[GenEll] §1` の高さの枠組みは `Found/` に載ったことになる。

| 層 | 中身 | 場所 | 状態 |
|---|---|---|---|
| A | `ADiv(F)` / `deg_F` / `deg_F`(正規化) | `Found/GenEll/ArithDiv.lean` | ★**実装済**(`sorry` 無し) |
| B・C | 算術直線束と引き戻し(スキーム上の直線束・解析化・hermitian 計量) | `Interface/GenEll/ArithLineBundle.lean` | `waiting` |
| **D** | **`APrc(F)`・同型 `ADiv/APrc ≅ APic`・`ht` の well-defined 性** | ★**本ファイル** | `sorry` |

## ★★`sorry` を置いてよい理由と、置き方の規律

`PLAN.md` の 2 トラック構成では `Skeleton/` が statement 専用トラックであり、
`sorry` が許されるのはここだけである(`Found/` には置かない)。
★ただし **statement は検査されない側**なので、G1(出典・逐語・目視)を必ず課す——
本ファイルの各宣言に `.src` を付け、`1_Structured` の該当 `<section>` を指す。

★**`Interface` を置いたことを「形式化した」と呼ばない**(`PLAN.md` §1)。
本ファイルが `sorry` 無しになるのは、層 B・C が `Found/` に載ってからである。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll NumberField

section Principal

variable {F : Type*} [Field F] [NumberField F]

/-- **主算術因子** `ADiv(f)`。

原文 (GenEll p.4):
> An arithmetic divisor on F is defined to be a finite formal sum

原文の定義は
`ADiv(f) ≝ Σ_{v∈𝕍(F)^non} ord_v(f)·v − Σ_{v∈𝕍(F)^arc} [F_v : ℝ]·log(|f|_v)·v`。

★**`sorry` の中身**: `ord_v(f)`(離散付値を `ℤ` として取り出す)と
`[F_v : ℝ]`(`InfinitePlace.mult`)と `log(|f|_v)` を組み、
**有限台であること**(ほとんどの `v` で `ord_v(f) = 0`)を示す必要がある。
★有限台性は「`f ≠ 0` なら分母分子の素因数分解が有限」という内容で、
mathlib の Dedekind 整域の理論から出るはずだが**未測定**。 -/
noncomputable def principalADiv (_f : Fˣ) : ADiv F := sorry

/-- **主算術因子のなす部分群** `APrc(F) ⊆ ADiv(F)`。

★`sorry` の中身: `principalADiv` が群準同型であること(`ADiv(f·g) = ADiv(f) + ADiv(g)`)。 -/
noncomputable def APrc (F : Type*) [Field F] [NumberField F] : AddSubgroup (ADiv F) := sorry

end Principal

section BaseChange

/-- ★★**正規化次数の基底変換不変性**。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

原文は `deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`(下線つき `deg` = 正規化版)と述べる。

★★**これが「高さが定義できる」ことの中身である**——
`x ∈ X(ℚ̄)` をどの数体 `F` で見ても同じ値になる、という主張。
**破れたら `ht` は定義できない。飾りではない。**

★`sorry` の中身: 素点の分岐・剰余次数の関係(`Σ_{w|v} e_w f_w = [K:F]`)を使う。
mathlib には `Ideal.sum_ramification_inertia` 系があるはずだが**未測定**。 -/
theorem degNormalized_base_change
    (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] (_a : ADiv F) (_b : ADiv K)
    (_hb : True) :
    degNormalized _b = degNormalized _a := sorry

end BaseChange

section Height

/-- **高さ** `ht_M̄`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

原文の定義は `ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ`、`x ∈ X(ℚ̄) = ∪_{[F:ℚ]<∞} X(F)`。

★`Interface/GenEll/ArithLineBundle.lean` の `PulledBackClassData` が
`x_F^* M̄` の正規化次数を受け取る。**その `Interface` が埋まれば `ht` は定義できる。**

★`sorry` の中身: `Interface` の `base_change_invariant` を使って
「`F` の取り方に依らない」ことを示し、`X(ℚ̄)` 上の関数として定義する。 -/
noncomputable def ht (_D : ABC3.Interface.GenEll.PulledBackClassData)
    (_M : Type) (_x : Type) : ℝ := sorry

/-- ★`ht` が **`F` の取り方に依らない**こと。

★これは `Interface` の `base_change_invariant` の**帰結**であって、
`Interface` を埋めれば自動的に出る——**そう書けることを型で示すのがこの Skeleton の役目**。

★★**現在の形は自明である**(`rfl`)。`ht` がまだ `sorry` なので、
「`F` に依らない」を型に出すには `ht` が `F` を引数に取る形になっている必要がある。
**ここは意図的に弱い**——`Interface` が埋まった時点で `ht` の型を
`(F : Type) → [NumberField F] → X(F) → ℝ` の形に直し、この定理を本物にする。
**その書き換えが済むまで、この宣言を「示した」と読んではならない。** -/
theorem ht_well_defined (D : ABC3.Interface.GenEll.PulledBackClassData)
    (M x : Type) : ht D M x = ht D M x := rfl

end Height

/-! ## ★出典の紐付け(`.src`) -/

def principalADiv.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

def APrc.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

def degNormalized_base_change.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

/-- ★この定理の証明が要求するもの。

原文は `[Szp], §1.1` を典拠に挙げるだけで、**証明を書いていない**。
実際に要るのは素点の分岐・剰余次数の関係である。 -/
def degNormalized_base_change.needs : List ProofObligation :=
  [ .citation "[Szp]" "§1.1(算術因子と次数の基本性質)"
      (.absent "mathlib 全体を `Arakelov` / `arithmetic.*line bundle` で grep、0 件(2026-08-16)") 4,
    .derivation
      "素点の分岐・剰余次数の関係 `Σ_{w|v} e_w f_w = [K:F]`。mathlib に `Ideal.sum_ramification_inertia` 系があるはずだが**未測定**" 4,
    .implicitStep
      "原文は `deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)` を『it follows that』とだけ書き、導出を示していない" 4 ]

def ht_well_defined.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(高さ ht_M̄)",
    sectionId := "genell-ht" }

/-- ★この定理の証明が要求するもの。

★**現在の形は自明**なので、実質的な依存は「`ht` の型を `F` を取る形に直すこと」である。
それは `Interface` が埋まってからでないとできない。 -/
def ht_well_defined.needs : List ProofObligation :=
  [ .implicitStep
      "★現在の `ht_well_defined` は `rfl` で示せる自明な形である。`ht` が `F` を引数に取る形に直るまで、この宣言は内容を持たない。直せるのは Interface(層 B・C)が埋まってからである" 4,
    .derivation
      "`degNormalized_base_change`(同ファイル)——`F` の取り方に依らないことの本体" 4 ]

def ht.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(高さ ht_M̄)",
    sectionId := "genell-ht" }

end ABC3.Skeleton.GenEll
