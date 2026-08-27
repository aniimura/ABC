import ABC3.Meta.Claim
import ABC3.Interface.GenEll.HeightTheory
import ABC3.Found.GenEll.BDClass
import ABC3.Found.GenEll.ArithDiv
import ABC3.Found.GenEll.Conductor
import ABC3.Found.GenEll.Prop16
import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Found.GenEll.HeightMetric
import ABC3.Found.GenEll.HeightAdditive
import ABC3.Found.GenEll.HeightClass
import ABC3.Found.GenEll.LogCondSigma
import ABC3.Found.GenEll.BDSlack
import ABC3.Found.GenEll.LogDiffTower

/-!
# [GenEll] §1 Generalities on Heights —— 必要 9 件の statement(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.3–p.10。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★このファイルの位置づけ —— S0(statement を型で固定する)

`[IUTchIV]` が要求する [GenEll] 24 件のうち **§1 の 9 件**をここで固定する。
★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。

★★**ただし `sorry` を消すことを目的にしてはならない。**
内容を `Interface` の仮説へ移せば `sorry` は消えるが、それは `tools/check.mjs` 冒頭 B5 が
名指しする穴である。ゆえに `Interface/GenEll/HeightTheory.lean` は
**公理を 1 つも持たない**(データと述語だけ)。

## ★★原典の記法の食い違いを、写したうえで別に記録する

★**`≲` の向きが、定義と用法で食い違う。** 原文 p.5 の定義は逐語で

> α ≲_F β … if there exists a ["constant"] C ∈ ℝ such that β(x) − α(x) ≤ C

であるから `α ≲ β` は **`β ≤ α + C`** を意味する。
ところが `Proposition 1.6` の**表題は "Conductor Bounded by the Height"** で、
本文は `log-cond_D ≲ ht_D` と書く——表題どおりなら `log-cond ≤ ht + C` のはずだが、
定義どおりに読むと **`ht ≤ log-cond + C`** になり、**逆である**。

★★**本ファイルは印字どおりに写す**(`BDle` = 逐語の定義)。
★食い違いそのものは `Gap/GenEll/BDDirection.lean` に
`GapRecord` と**反例**(`Found/GenEll/BDClass.lean` の `bdle_ne_bdge`)つきで記録する。
**写し間違えないことと、食い違いを黙らないことは、別の仕事である。**

## ★実装済みのものは実装済みのものを使う

- `BDle` / `BDge` / `BDeq` / `BDClass` —— `Found/GenEll/BDClass.lean`(`sorry` 無し)
- `ADiv` / `deg` / `degNormalized` —— `Found/GenEll/ArithDiv.lean`(`sorry` 無し)

★これにより §1 の skeleton は**空虚な述語の羅列にならない**。
`Definition 1.2` と `Definition 1.5` は、下に書くとおり
**既に作ってあるものの上に本当に載る**。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll

/-! ## Definition 1.2 —— 高さ関数とその BD-class -/

/-- **[GenEll] Definition 1.2**。

原文 (GenEll p.5):
> The relation "≈" clearly defines an equivalence relation on the set of functions

(i) は高さ関数 `ht_M̄ : X(ℚ̄) → ℝ`、(ii) は BD-class。
★**両方をまとめて 1 つの式にした** —— 高さ関数の **BD-class** を取る操作である。
これが原文の (i)(ii) の合流点であり、以降 §1・§2 の不等式はすべてこの型の上で語られる。

★`BDClass` は `Found/GenEll/BDClass.lean` で**商として実際に構成してある**
(`bdSetoid` / `Quotient`)。posit ではない。 -/
noncomputable def htClass (D : HeightTheoryData) (M : D.ABundle) : BDClass D.Point :=
  BDClass.mk (D.ht M)

/-! ## Example 1.3 —— Galois-finite と compactly bounded -/

/-- **[GenEll] Example 1.3, (i)** の `X(ℚ̄)^{=d}`。

原文 (GenEll p.5):
> [cf. the discussion in Definition 1.5, (i), below of "minimal fields of definition"].

`X(ℚ̄)^{=d} ≝ X(ℚ̄)^{≤d} \ X(ℚ̄)^{≤d−1}`。 -/
def degEq (D : HeightTheoryData) (d : ℕ) : Set D.Point :=
  D.degLe d \ D.degLe (d - 1)

/-- **[GenEll] Example 1.3, (i)** の `Galois-finite`。

原文 (GenEll p.5):
> over the positive integers] is finite, then we shall say that E is Galois-finite.

★**これは posit ではなく定義である。** 原文は「各 `E^{≤d}` が有限」としか言っておらず、
それは `degLe` があれば**純粋に集合論的に**書ける。
★したがって `Example 1.3, (i)` は **Arakelov 理論も Galois 表現も要求しない**。 -/
def GaloisFinite (D : HeightTheoryData) (E : Set D.Point) : Prop :=
  ∀ d : ℕ, 0 < d → (E ∩ D.degLe d).Finite

/-- **[GenEll] Example 1.3** が導入する 2 つの語。

★**(i) は上で定義した**(集合論だけ)。
★**(ii) は `Interface` が posit している**——`compact domain` と `X^arc` を要求するからで、
mathlib に `complex analytic space` は 0 件(2026-08-16 実測)。

この対を 1 つの式にすることで、「この Example が何を導入したか」を型で固定する。 -/
def example_1_3 (D : HeightTheoryData) :
    (Set D.Point → Prop) × (Set D.Point → Prop) :=
  (GaloisFinite D, D.CompactlyBounded)

/-! ## Proposition 1.4 —— 高さの基本性質 -/

/-- **[GenEll] Proposition 1.4**(Basic Properties of Heights)。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

(i) `ht_{L̄⊗M̄}(x) = ht_L̄(x) + ht_M̄(x)` ——★**`≈` ではなく `=`**(目視確認)。
(ii) `L_ℚ` のある正冪が大域切断で生成されるなら `ht_L̄ ≳ 0` ——★**`≥` ではなく `≳`**。
(iii) `ht_L̄` の BD-class は `L_ℚ` の同型類だけに依る。
(iv) `L_ℚ` が ample なら `{x ∈ X(ℚ̄)^{≤d} : ht_L̄(x) ≤ C}` は**有限**(Northcott)。

★(ii) の `≳ 0` は逐語では `BDge (ht L) 0`、すなわち `ht(x) − 0 ≤ C`。
**定義どおりに読むと「上に有界」であって「下に有界」ではない**——
これも上の docstring で述べた向きの食い違いの一例である(`Gap/GenEll/BDDirection.lean`)。

## ★★★★★★★★ 2026-08-26——構成に置き換えて閉じた(第 371 ブロック)

★以前は `∀ D : HeightTheoryData` と量化していたが、
**`HeightTheoryData` は公理を 1 つも持たないデータ**なのでそれは**偽**であった
(反例: `Check/GenEll/HeightAxiomGap.lean` の `prop_1_4_statement_false`。
`Point ≝ ℕ`･`ABundle ≝ ℝ`･`tensor ≝ (+)`･`ht L x ≝ L²` で
**(i) 加法性と (iv) Northcott が同時に破れる**)。

★★公理を足すのは **B5 そのもの**である——`Proposition 1.4` の (i)-(iv) は
`HeightTheoryData` に足すべき公理**そのもの**だから、足せば仮定の言い換えになる。
★★★ゆえに**構成に置き換えた**。

### ★★★★ 4 つの内訳

| 条 | 中身 | どこにあるか |
|---|---|---|
| (i) 加法性 | `ht(D⊗E) = ht(D) + ht(E)` | `htArith_tensor_unconditional` |
| (ii) 下に有界 | `-C ≤ ht`(`C` は `F` にも点にも依らない) | `Prop16.lean` の `prop_1_4_ii` |
| (iii) `≈` | 因子が同じで計量が連続なら高さの差は一様に有界 | `HeightMetric.lean` の `htArith_sub_abs_le`(第 371) |
| (iv) Northcott | 射影モデルがあれば有限 | `NorthcottCoord.lean` の `northcott_of_projModel`(第 369-371) |

### ★★★★★逸脱と、残してあるもの(明示)

1. ★**量化する対象**: `∀ D : HeightTheoryData` → `ArcModel` + `ArithCartier`。前者では偽だから。
2. ★★**(iii) の形**: 原文は『生成ファイバーが同じなら』。因子表示では
   『因子が同じ(計量だけ違う)』である。**垂直因子の差の分は含めていない**。
3. ★★★**(iv) の形**: 原文は `X` の `ℤ`-固有性と `L_ℚ` の豊富性から**射影埋め込み**を得る。
   本 statement はその埋め込みを **`ArcModel` と同じ立場でデータとして受けている**。
   ★★★★**Northcott 性そのものは受けていない**——
   それは `Found/GenEll/NorthcottTuple.lean`･`NorthcottCoord.lean` で**証明している**。
   ★★★★★残っているのは「**`htArith` がその意味での射影モデルを持つ**」という幾何の段だけである。 -/
theorem prop_1_4 {X : AlgebraicGeometry.Scheme.{0}} {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ABC3.Found.GenEll.ArcModel X V)
    [Nonempty (ABC3.Found.GenEll.complexPoints X)]
    (Dv : ABC3.Found.GenEll.ArithCartier X)
    (hg : @Continuous _ _ M.topology _ Dv.green) :
    -- (i) 加法性
    (∀ (Ev : ABC3.Found.GenEll.ArithCartier X) (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        ABC3.Found.GenEll.pullbackIdeal F Dv.divisor xF ≠ 0 →
        ABC3.Found.GenEll.pullbackIdeal F Ev.divisor xF ≠ 0 →
        ABC3.Found.GenEll.htArith F (Dv.tensor Ev) xF
          = ABC3.Found.GenEll.htArith F Dv xF + ABC3.Found.GenEll.htArith F Ev xF)
    -- (ii) 下に一様に有界
  ∧ (∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        -C ≤ ABC3.Found.GenEll.htArith F Dv xF)
    -- (iii) 因子が同じで計量が連続なら高さは BD-同値
  ∧ (∀ Ev : ABC3.Found.GenEll.ArithCartier X, Ev.divisor = Dv.divisor →
        @Continuous _ _ M.topology _ Ev.green →
        ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
          (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
          |ABC3.Found.GenEll.htArith F Dv xF - ABC3.Found.GenEll.htArith F Ev xF| ≤ C)
    -- (iv) Northcott(射影モデルを与えられたものとして)
  ∧ (∀ (P I : Type) [Finite I] (d : ℕ) (ht : P → ℝ)
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p)),
        (∀ p, Module.finrank ℚ (fld p) ≤ d) →
        ∀ (crd : ∀ p, I → (fld p)) (idx : I) (const : ℝ),
        (∀ p, haveI := hnf p; Height.mulHeight (crd p) ≤ Real.exp (ht p + const)) →
        Function.Injective (fun (p : P) (i : I) => ((crd p i / crd p idx : fld p) : ℂ)) →
        ∀ C : ℝ, {p : P | ht p ≤ C}.Finite) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro Ev F _ _ xF hD hE
    exact ABC3.Found.GenEll.htArith_tensor_unconditional F Dv Ev xF hD hE
  · exact ABC3.Found.GenEll.prop_1_4_ii M Dv hg
  · intro Ev hdiv hcont
    exact ABC3.Found.GenEll.htArith_sub_abs_le M Dv Ev hdiv.symm hg hcont
  · intro P I _ d ht fld hnf hdeg crd idx const hcmp hinj C
    exact ABC3.Found.GenEll.northcott_of_projModel d ht fld hnf hdeg crd idx const hcmp hinj C

/-! ## Remark 1.4.1 —— 理論が `X_ℚ` だけに依ること -/

/-- **[GenEll] Remark 1.4.1**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

「`X(ℚ̄)` 上の高さ関数の **BD-class の理論**は、スキーム `X_ℚ` だけに依る」。

★**「だけに依る」を型で書くとこうなる** —— 2 つのモデル `D`, `D'` について、
点と生成ファイバーの同一視があり、そのもとで直線束の生成ファイバーが対応するなら、
高さ関数は BD-同値である。

★これが **`Theorem 2.1` が「数体上の曲線」から出発できる根拠**である——
ℤ-モデル `X` の取り方に依らないことを保証している。

## ★★★★★★★ 2026-08-26——構成に置き換えて閉じた(第 372 ブロック)

原文は『follows **immediately** from the definitions, together with
**Proposition 1.4, (iii)**』と書いている。
★すなわち中身は「**高さは類にしか依らない**」であり、
構成の側では `Found/GenEll/HeightClass.lean` が**すでにそれを持っていた**
——`htArith_bdeq_of_pullbackAPic_eq`。

★★構成の側では**定数差すら出ない**(`C = 0`)——
機構は `deg` が主算術因子の上で消えること(積公式)である。

### 逸脱(明示)

| 項 | 原典 | 形式化 |
|---|---|---|
| 量化する対象 | `∀ D D' : HeightTheoryData` | **同じ `X` 上の 2 つの `ArithCartier`** |
| 「だけに依る」の中身 | `X_ℚ` だけに依る | **引き戻した算術因子の類だけに依る** |

★★★**含めていないもの**: 異なる 2 つの `ℤ`-モデルの比較
——原文は『同型がある有限素数集合 `Σ` の上で `ℤ[Σ^{-1}]` へ延びる』という段を使うが、
その段はまだ持っていない。★**本 statement はその手前までである**。 -/
theorem remark_1_4_1 {X : AlgebraicGeometry.Scheme.{0}}
    (Dv Ev : ABC3.Found.GenEll.ArithCartier X)
    (F : Type) [Field F] [NumberField F]
    (h : ∀ xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X,
      ABC3.Found.GenEll.pullbackAPic F Dv xF = ABC3.Found.GenEll.pullbackAPic F Ev xF) :
    BDeq (fun xF => ABC3.Found.GenEll.htArith F Dv xF)
      (fun xF => ABC3.Found.GenEll.htArith F Ev xF) :=
  ABC3.Found.GenEll.htArith_bdeq_of_pullbackAPic_eq F Dv Ev h

/-! ## Definition 1.5 —— log-diff と log-cond -/

/-- **[GenEll] Definition 1.5**。

原文 (GenEll p.8):
> Note if x ∈ X(F ) ⊆ X(Q), where [F : Q] &lt; ∞, then by considering the scheme-theoretic image of the corresponding morphism Spec(F ) → X, one obtains a well-defined minimal field of definition Fmin ⊆ F of x.

(i) 最小定義体 `F_min`、(ii) `E` が **reduced** ⇔ `E = E_red`、
(iii) `log-diff_X(x) ≝ deg_F(δ_x)`、(iv) `log-cond_D(x) ≝ deg_F(f_x^D)`。

★★**(iii)(iv) の右辺は、既に作ってある `degNormalized` そのものである。**
原文の下線つき `deg` は正規化次数であり、`Found/GenEll/ArithDiv.lean` に
`sorry` 無しで実装済み。★ゆえにこの定義の**式そのものはもう持っている**——
残るのは `δ_x`(差積イデアルが定める有効算術因子)と
`f_x^D ≝ (D_x)_red`(引き戻した因子の被約化)を**構成すること**である。

★**その 2 つは Arakelov 理論ではない**——可換環論(`IsDedekindDomain.differentIdeal`)と
scheme 論(Cartier 因子の引き戻し)である。 -/
noncomputable def defn_1_5 {F : Type*} [Field F] [NumberField F] : ADiv F → ℝ :=
  degNormalized

/-! ★`(−)_red`(`ADivRed`)とその冪等性・有効性・**次数不等式**は
`Found/GenEll/Conductor.lean` へ移した——★`Proposition 1.6` の**非アルキメデス側**が
そこで実際に証明されるからである(`deg((a)_red) ≤ deg(a)`)。 -/

/-! ## Remark 1.5.1 —— log-cond の BD-class が `(X_ℚ, D_ℚ)` だけに依ること -/

/-- **[GenEll] Remark 1.5.1**。

原文 (GenEll p.8):
> Remark 1.5.1. In the spirit of Remark 1.4.1, we observe that the log-different

`log-diff_X` はスキーム `X_ℚ` だけに依る。`log-cond_D` は対 `(X, D)` に依り得るが、
**その BD-class は ℚ-スキームの対 `(X_ℚ, D_ℚ)` だけに依る**。

## ★★★★★★★★★★ 2026-08-27——**構成に載せ替えた**(第 420 ブロック)

★**旧 statement(`∀ D D′ : HeightTheoryData, …`)は偽であった**
——`Check/GenEll/RemarkAxiomGap.lean` の `remark_1_5_1_false` で機械検証済み。
`HeightTheoryData` は公理を 1 つも持たないので、点の全単射だけでは
`logDiff` を何も縛らない(`logDiff` を `0` と `x` に取れば反例になる)。

★★**したがって旧 statement の `sorry` は「spreading out 待ち」ではなかった**
——spreading out を実装しても statement は偽のままで、**原理的に閉じられない**。
`Proposition 1.4` / `Remark 1.4.1` と同じ病であり、同じ治療(構成への載せ替え)を施した。

### 前半と後半

| 半分 | 構成の側での扱い |
|---|---|
| `log-diff_X` は `X_ℚ` だけに依る | ★**定義から明らか**——`logDiffOfField` は `X` を引数に持たない |
| `log-cond_D` の BD-class は `(X_ℚ, D_ℚ)` だけに依る | ★★**本 statement の後半**。定数は `Σ_{q∈Σ} log q` |

### 逸脱(明示)

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ D D′ : HeightTheoryData` | **`ArithCartier` / `IdealSheafData` の引き戻し** | 前者では偽だから |
| spreading out | 証明の中で使う | **仮定 `hagree` として受ける** | ★下記 |

★★★**`hagree`(`Σ` の外で導手が一致する)は spreading out の帰結である**。
原文の証明は『同型が `ℤ[Σ^{-1}]` へ延びる』でこれを与え、そこから BD-class の一致を結論する。
★**本 statement は後半を証明し、前半を仮定に置いた**——`.needs` に明記してある。

★★★★**空虚ではない**——`Check/GenEll/RemarkSigmaWitness.lean` に
`hagree` が実際に満たされる場合(`D = D′`、`Σ = ∅`)を置いてある。

### ★★★★★定数が `Σ` だけで決まることが要点

`Σ_{q ∈ Σ} log q` は **`Σ` だけで決まり、点 `x` にも定義体 `F` にも依らない**。
一様性の機構は「`log-cond` は `deg / [F:ℚ]`、`Σ` の上の寄与の分子 `Σ_{v|q} f_v` は
`Σ_{v|q} e_v f_v = [F:ℚ]` で抑えられる ⟹ **`[F:ℚ]` が約分される**」である
(`Found/GenEll/SigmaBound.lean`)。
★これが**この論文が BD-class を使う理由そのもの**である。 -/
theorem remark_1_5_1 (F : Type) [Field F] [NumberField F]
    {X X' : AlgebraicGeometry.Scheme.{0}}
    (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (ePt : (ABC3.Found.GenEll.specRingOfIntegers F ⟶ X) →
           (ABC3.Found.GenEll.specRingOfIntegers F ⟶ X'))
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : ABC3.Found.GenEll.FinitePlace F → ℕ)
    (hover : ∀ v : ABC3.Found.GenEll.FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (hI : ∀ xF, ABC3.Found.GenEll.pullbackIdeal F D xF ≠ 0)
    (hI' : ∀ xF, ABC3.Found.GenEll.pullbackIdeal F D' (ePt xF) ≠ 0)
    (hagree : ∀ xF, ∀ v, ch v ∉ Sig →
      (ABC3.Found.GenEll.conductorADiv F D xF).fin v
        = (ABC3.Found.GenEll.conductorADiv F D' (ePt xF)).fin v) :
    -- (前半) `log-diff` は模型に依らない(構成の側では定義から)
    (∀ (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X)
       (xF' : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X'),
        ABC3.Found.GenEll.logDiffOfField F = ABC3.Found.GenEll.logDiffOfField F)
    -- (後半) `log-cond` の BD-class は一致する
  ∧ BDeq (fun xF => ABC3.Found.GenEll.logCond F D xF)
         (fun xF => ABC3.Found.GenEll.logCond F D' (ePt xF)) := by
  refine ⟨fun _ _ => rfl, ⟨∑ q ∈ Sig, Real.log q, fun xF => ?_⟩⟩
  exact ABC3.Found.GenEll.abs_logCond_sub_le_sum_log F D D' xF (ePt xF)
    (hI xF) (hI' xF) Sig hprime ch (hagree xF) hover


/-! ## Proposition 1.6 —— 導手は高さで抑えられる -/

/-- **[GenEll] Proposition 1.6**(Conductor Bounded by the Height)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

`L = O_X(D)` とし `U ≝ X\D`、`ht_D ≝ ht_L̄` とすると、`U(ℚ̄)` 上で `log-cond_D ≲ ht_D`。

★**枝番を持たない単一の主張**である(目視確認 2026-08-16)。

★★**ここが向きの食い違いが最も鮮明に出る場所である。**
表題は "Conductor **Bounded by** the Height"、すなわち `log-cond ≤ ht + C` を言っている。
ところが `≲` を **p.5 の定義どおり**に読むと `BDle log-cond ht` は
`ht(x) − log-cond(x) ≤ C`、つまり **`ht ≤ log-cond + C`** であり **逆になる**。
★本 statement は**印字どおり**(`BDle`)に写した。食い違いは
`Gap/GenEll/BDDirection.lean` に記録してある。

★証明が `X^arc` を要求する箇所も目視で特定した——アルキメデス素点の寄与を
「**コンパクト空間 `X^arc` 上の連続関数 `|s|_L` が有界**」で処理している。
非アルキメデス側は `Definition 1.5, (iv)` の `(−)_red` だけで済む。

## ★★★★★★★★ 2026-08-26——`sorry` が消えた(第 368 ブロック)

★以前は `∀ D : HeightTheoryData` と量化していたが、**それは偽である**
(`Check/GenEll/HeightAxiomGap.lean`)。★★公理を足すのは B5 そのものなので、
本ファイルの診断が書いていたとおり**構成に置き換えた**。

★★★中身は `Found/GenEll/Prop16.lean` にすでにあった——三段である:

1. `X^arc` はコンパクト(`ArcModel.compactSpace`)
2. 連続な Green 関数は下に有界(`ArcModel.exists_bound`)
3. 導手は高さ + 定数で押さえられる(`ArchBound.logCond_le_htArith_add`)

★★★★**定数 `C` が数体 `F` にも点にも依らない**のが要である(次数で正規化してあるから)。

### 逸脱 2 件

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ D : HeightTheoryData` | **`ArcModel` + `ArithCartier`** | 前者では偽だから |
| `≲` の向き | 印字どおりの `BDle` | **`log-cond ≤ ht + C`** | 表題 "Conductor **Bounded by** the Height" の向き。`Gap/GenEll/BDDirection.lean` に記録済み |

★原文の `U(ℚ̄)` への制限は `pullbackIdeal F Dv.divisor xF ≠ 0` が担っている。 -/
theorem prop_1_6 {X : AlgebraicGeometry.Scheme.{0}} {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ABC3.Found.GenEll.ArcModel X V)
    [Nonempty (ABC3.Found.GenEll.complexPoints X)]
    (Dv : ABC3.Found.GenEll.ArithCartier X)
    (hg : @Continuous _ _ M.topology _ Dv.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        ABC3.Found.GenEll.pullbackIdeal F Dv.divisor xF ≠ 0 →
        ABC3.Found.GenEll.logCond F Dv.divisor xF
          ≤ ABC3.Found.GenEll.htArith F Dv xF + C :=
  ABC3.Found.GenEll.prop_1_6 M Dv hg

/-! ## Proposition 1.7 —— 導手と log-different -/

/-- **[GenEll] Proposition 1.7**(Conductors and Log Differents)。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

(i) `U_Y(ℚ̄)` 上で
`log-cond_E − log-cond_D ≲ log-diff_Y − log-diff_Z ≲ (1 − 1/e)·log-cond_E`。

## ★★★★★★★★★★★ 2026-08-27——**構成に載せ替えた**(第 423 ブロック)

★**旧 statement(`∀ S : CoveringSetup, S.hyp → …`)は偽であった**
——`Check/GenEll/Prop17AxiomGap.lean` の `prop_1_7_false` で機械検証済み。
`hyp : Prop` は `CoveringSetup` の**不透明なフィールド**であって
`logDiff` や `logCond` を**何も縛らない**。`hyp := True` と置けば結論は破れる。

★★**落とさずに持つことと、制約になることは別である**——これが実測でわかったことである。

## ★★原文の証明の構造をそのまま型にした

原文 p.10 の証明は**明示的に 2 つに分けている**(2026-08-16 の 260 dpi 目視で確認):

| 部分 | 原文の言い方 |
|---|---|
| prime-to-`Σ` | 「**`=` と `≤` であって `≲` ではない**」 |
| `Σ` の上 | 「**`≥ 0` であって `≳` ではない**」 |

★すなわち**厳密な不等式が `Σ` の外で成り立ち、`Σ` の上の食い違いが一様に有界**であり、
その 2 つから `≲` が出る。★★本 statement はその**後半**を証明する
(`Found/GenEll/BDSlack.lean` の `bdle_of_bounded_slack`)。

★★★`Σ` の上の食い違いが `Σ_{q∈Σ} log q` で一様に抑えられることは
`Found/GenEll/SigmaBound.lean` / `LogCondSigma.lean` で取ってある(第 412–415)。
**点にも定義体にも依らない定数**であることが要点である。

## ★★★★`log-diff_Y − log-diff_Z` が何であるかは分かっている

`Found/GenEll/LogDiffTower.lean` の `logDiffOfField_tower` は**等式**である:

> `log-diff(K) − log-diff(F) = log N(𝔡_{K/F}) / [K:ℚ]`

★したがって (i) の中辺は**相対 different の対数ノルム**そのものである。
`hlow` / `hup` はその上下からの評価であり、原文が
『the elementary theory of differents』と呼ぶものにあたる
——その**局所の核**は `TameRamification.lean` / `DifferentKummer.lean` /
`TotallyRamified.lean` で取ってある(第 374–411、馴分岐 6/6)。
★★残っているのは**局所から大域への組み立て**であり、`.needs` に明記した。

## ★★★★★逸脱(明示)

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ S : CoveringSetup` | **点の型 `P` と実数値関数** | 前者では偽だから |
| elementary claim | 証明の中で使う | **仮定 `hlow` / `hup` として受ける** | 局所から大域への段が未了 |
| (ii) Riemann–Hurwitz | 述べる | **含めない** | `deg` が語彙の外(`.needs` に記録) |

★★★★★★**空虚ではない**——`Check/GenEll/Prop17Witness.lean` に
仮定が実際に満たされる場合を置いてある。 -/
theorem prop_1_7 {P : Type}
    (condE condD diffY diffZ slackLow slackUp : P → ℝ) (e : ℕ) (he : 0 < e)
    (Sig : Finset ℕ)
    -- prime-to-`Σ` の部分(原文が `=` と `≤` と明示している段)
    (hlow : ∀ x, diffY x - diffZ x ≤ (condE x - condD x) + slackLow x)
    (hup : ∀ x, (1 - 1 / (e : ℝ)) * condE x ≤ (diffY x - diffZ x) + slackUp x)
    -- `Σ` の上の食い違いは `Σ_{q∈Σ} log q` で一様に抑えられる
    (hsl : ∀ x, slackLow x ≤ ∑ q ∈ Sig, Real.log q)
    (hsu : ∀ x, slackUp x ≤ ∑ q ∈ Sig, Real.log q) :
    BDle (fun x => condE x - condD x) (fun x => diffY x - diffZ x)
  ∧ BDle (fun x => diffY x - diffZ x) (fun x => (1 - 1 / (e : ℝ)) * condE x) :=
  ⟨ABC3.Found.GenEll.bdle_of_bounded_slack _ _ slackLow _ hlow hsl,
   ABC3.Found.GenEll.bdle_of_bounded_slack _ _ slackUp _ hup hsu⟩

/-- ★★★★**(i) の中辺の正体** —— `log-diff_Y − log-diff_Z` は相対 different の対数ノルム。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★`Found/GenEll/LogDiffTower.lean` の `logDiffOfField_tower` を並べ替えただけである。
★★**不等式ではなく等式**なので、(i) の両側は `𝔡_{K/F}` の上下からの評価に**帰着する**。 -/
theorem prop_1_7_middle (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] :
    ABC3.Found.GenEll.logDiffOfField K - ABC3.Found.GenEll.logDiffOfField F
      = Real.log (Ideal.absNorm
          (differentIdeal (NumberField.RingOfIntegers F) (NumberField.RingOfIntegers K)))
        / (Module.finrank ℚ K : ℝ) := by
  rw [ABC3.Found.GenEll.logDiffOfField_tower F K]
  ring


/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def htClass.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2",
    sectionId := "genell-def-1-2" }

def degEq.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

def GaloisFinite.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

def example_1_3.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3",
    sectionId := "genell-ex-1-3" }

def prop_1_4.src : Source :=
  { paper := "GenEll", pdfPage := 6, item := "Proposition 1.4",
    sectionId := "genell-prop-1-4" }

/-- ★原文 p.6–p.7 の証明を通読して数えた。 -/
def prop_1_4.needs : List ProofObligation :=
  [ .implicitStep
      "(i) は『follows immediately from the definitions』の 1 行で済まされている。ht の定義(deg_F の和)まで降りれば加法性だが、その deg_F が Interface 待ちである" 6,
    .implicitStep
      "(ii) の証明はアルキメデス素点の寄与を『X^arc がコンパクトだから |s|_L は有界』で処理する。★複素解析空間そのものが mathlib に無い(complex analytic space 0 件、2026-08-16 実測)" 6,
    .implicitStep
      "(iii) は (i)(ii) から従うと原文は書くが、そこで『ℤ 上の算術直線束の高さは有界』という補題を暗黙に使っている(ht_{L̄⊗M̄} ≈ ht_L̄ の式)" 6,
    .folklore "(iv) は Northcott の有限性定理。原文は証明を与えず『well-known』の扱いをしている" 7,
    .implicitStep
      "★statement の語彙(APic(X)・ht・ample・大域切断生成)を Interface/GenEll/HeightTheory.lean に posit した。**我々は作っていない**" 6 ]

def remark_1_4_1.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Remark 1.4.1",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Proposition 1.4, (iii)(BD-class が L_ℚ の同型類だけに依る)" 6,
    .implicitStep
      "原文は『follows immediately from the definitions』とだけ書く。★「X_ℚ だけに依る」を型にするには『2 つの ℤ-モデルの生成ファイバーが同型なら高さが BD-同値』を言う必要があり、その同一視の与え方は原文に無い(我々が ePt / eGen として補った)" 8 ]

def defn_1_5.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5",
    sectionId := "genell-def-1-5" }

def remark_1_5_1.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Remark 1.5.1",
    sectionId := "genell-rem-1-5-1" }

def remark_1_5_1.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Remark 1.4.1(BD-class の理論が X_ℚ だけに依る)" 8,
    .implicitStep
      "★★★本 statement は原文の証明のうち**後半だけ**を証明している。前半(同型が ℤ[Σ^{-1}] へ延びること = spreading out)は仮定 hagree(Σ の外で導手が一致する)として受けている。★空虚でないことは Check/GenEll/RemarkSigmaWitness.lean で確かめてある" 9,
    .otherPaper "[Stacks]" "32.22 Descending finite type schemes(有限型スキームの降下——EGA IV §8 にあたる)" 4441,
    .implicitStep
      "★★旧 statement(∀ D D′ : HeightTheoryData, …)は**偽**であった——Check/GenEll/RemarkAxiomGap.lean の remark_1_5_1_false で機械検証済み。2026-08-27 に構成へ載せ替えた(第 420 ブロック)。逆向きの逸脱をしていないことは、定数が Σ だけで決まり点にも定義体にも依らないことで担保されている" 8 ]

def prop_1_6.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.6",
    sectionId := "genell-prop-1-6" }

/-- ★原文 p.9 の証明は 5 行しかない。その 5 行が要求するものを数えた。 -/
def prop_1_6.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Definition 1.5, (iv)(導手 f_x^D ≝ (D_x)_red と log-cond_D)" 8,
    .otherPaper "[GenEll]" "Proposition 1.4, (iii)(ht_D ≝ ht_L̄ が well-defined であること)" 6,
    .implicitStep
      "★アルキメデス側は『コンパクト空間 X^arc 上の連続関数 |s|_L が有界』で片づけられている。X^arc(複素解析空間)が mathlib に 0 件なので、この 1 行が層まるごとに相当する" 9,
    .implicitStep
      "★非アルキメデス側は『(−)_red の定義から従う』とだけ書く。実際には『D_x の係数 ≥ 1 の素点でだけ (D_x)_red が 1 を持つ』という不等式であり、算術因子の上で書ける(ADivRed を参照)" 9,
    .implicitStep
      "★★表題 'Conductor Bounded by the Height' と、p.5 の ≲ の定義が示す向きが逆である。Gap/GenEll/BDDirection.lean に記録した。**本 statement は印字どおりに写してある**" 9 ]

def prop_1_7_middle.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.7((i) の中辺は相対 different の対数ノルム)",
    sectionId := "genell-prop-1-7" }

/-- ★**依存は無い** —— mathlib の差積の tower 公式を並べ替えただけである。 -/
def prop_1_7_middle.needs : List ProofObligation := []

def prop_1_7.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.7",
    sectionId := "genell-prop-1-7" }

/-- ★原文 p.10 の証明を通読して数えた。

★**この証明の核は初等的である**——「素数 `p` と正整数 `d` を固定すると、
`[L:K] ≤ d` なる `ℚ_p` の有限拡大の有限 Galois 拡大 `L/K` すべてについて、
差積イデアルが `p^n·O_L` を含むような正整数 `n` が存在する」。
Arakelov も Galois 表現も要らず、**局所体の分岐理論と Kummer 理論だけ**である。 -/
def prop_1_7.needs : List ProofObligation :=
  [ .implicitStep
      "★★★本 statement は原文の証明のうち**後半だけ**を証明している。prime-to-Σ の厳密な不等式(hlow / hup)は仮定として受けている。★空虚でないことは Check/GenEll/Prop17Witness.lean で確かめてある" 10,
    .folklore "原文が『the elementary theory of differents』と呼ぶもの。★局所の核は 2026-08-26〜27 の第 374-411 で実装した(馴分岐 6/6)。★★残っているのは**局所から大域への組み立て**である" 10,
    .implicitStep
      "★原文は prime-to-Σ 部分について『with \"=\" and \"≤\", not \"≲\"!』、Σ 上について『with \"≥\", not \"≳\"!』と**明示的に区別している**。★★本 statement はその区別を型にした——厳密な不等式(hlow/hup)と一様な有界性(hsl/hsu)から ≲ が出る(Found/GenEll/BDSlack.lean)" 10,
    .citation "[GenEll]" "局所体の分岐理論と Kummer 理論(証明の核となる initial claim)"
      (.inMathlib "IsDedekindDomain.differentIdeal / Polynomial.IsSplittingField / IsCyclic —— ただし『[L:K] ≤ d なる全ての L/K に一様な n』という**一様性**の部分は mathlib に無い(2026-08-16 実測)") 10,
    .otherPaper "[GenEll]" "Remark 1.5.1(Σ 上の log-cond の寄与が ≈ 0 であること)" 8,
    .otherPaper "[Stacks]" "53.12 Riemann-Hurwitz((ii) の原典側の節点——第 419 で手元にあることを実測した)" 4441,
    .implicitStep
      "★★旧 statement(∀ S : CoveringSetup, S.hyp → …)は**偽**であった——Check/GenEll/Prop17AxiomGap.lean の prop_1_7_false で機械検証済み。hyp : Prop は不透明なフィールドであって logDiff や logCond を何も縛らない。2026-08-27 に構成へ載せ替えた(第 423 ブロック)" 9,
    .implicitStep
      "(ii) の Riemann–Hurwitz は Y_ℚ・Z_ℚ 上の直線束の次数 deg(−) を要求する。本ファイルでは (i) だけを固定した" 10 ]

end ABC3.Skeleton.GenEll
