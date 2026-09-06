import ABC3.Found.PGC.LubinTateEndoLimit
import ABC3.Found.PGC.LubinTateClosure
import ABC3.Found.PGC.LubinTateActionInclusion

/-!
# `K_π` は `f` に依らない —— 同じ `π` の `f, g ∈ F_π` で `K_{π,f} = K_{π,g}`(`sorry` 無し)

経路 Λ の節点 Λ6a′。`Found/PGC/LubinTateClosure.lean` が「逸脱(記録)」
として明示的に**主張しないまま残していた**独立性

> `f` への依存: `lubinTateClosure` は一意化元 `π` だけでなく Lubin-Tate
> 級数 `f` にも依存する形で定義されている。古典論では `K_π` は `f` に
> 依らないが、その独立性は Λ6(Lubin-Tate-Rosen)の担当なので、ここでは
> 主張しない。

の**うち `π` を固定した部分**を閉じる。すなわち、同じ素元 `π` に属する
2 つの Lubin-Tate 級数 `f, g ∈ F_π` に対し

`lubinTateClosure(f) = lubinTateClosure(g)`

を示す。★`π` そのものを取り替える場合(`π ≠ π'`)は別の塔になり得て、
そちらは Dwork の補題(`Found/PGC/DworkTheta.lean` 系)を通る本体
Λ6 の担当である——本ファイルは Dwork を一切通らない。

## 典拠

Milne, `Class Field Theory`, Prop. 3.10 の **`u = 1` の場合**。原典は
一般の単数 `u ∈ 𝒪_K^×` について `[u]_{g,f} ∈ 𝒪_{K̂^ur}[[T]]` を作り、
`σθ = θ ∘ [u]` を満たす `θ` を Dwork の補題で `𝒪_K` 係数へ降ろす。
`u = 1` なら `ε = 1` が取れて `σθ = θ`、すなわち `θ` は最初から
`𝒪_K[[T]]` の元 `[1]_{g,f}` であり、`K̂^ur`・`ℂ_K`・Dwork のどれも要らない。

★`.src` は書いていない。典拠の Milne CFT は `ResearchPaper/1_Structured/`
に構造化されていないため、`ABC3.Meta.Source` の必須フィールド `sectionId`
を正直に埋められない(嘘の `sectionId` は書かない)。

## 何を足したか

### 抽象核(Lubin-Tate・分岐・Galois の語彙が 1 つも出てこない)

一般の可換環 `A` 上の形式冪級数だけの主張:

| 宣言 | 内容 |
|---|---|
| `constantCoeff_iteratedLubinTate_of_constantCoeff_zero` | 定数項 `0` の級数の `n` 回自己合成も定数項 `0` |
| `hasSubst_iteratedLubinTate` | 同上を `HasSubst` の形で |
| `subst_iteratedLubinTate_of_intertwine` | `f ∘ θ = θ ∘ g` ならば `f^{(n)} ∘ θ = θ ∘ g^{(n)}`(合成の結合律だけ) |
| `subst_eq_X_of_intertwine_pair` | `f∘θ=θ∘g`・`g∘η=η∘f`・`θ'(0)=η'(0)=1` ならば `θ ∘ η = X` |

★`subst_iteratedLubinTate_of_intertwine` は `CommRing A` だけで成り立つ
(局所環すら要らない)。`subst_eq_X_of_intertwine_pair` だけが一意性補題
`powerSeries_uniqueness` を経由するので `IsLocalRing`・`IsDomain`・
`maximalIdeal = span {π}` を要求する——それでも付値・分岐の語彙は出ない。

### 具体層

| 宣言 | 内容 |
|---|---|
| `norm_lubinTateEvalAtPoint_le` | 定数項 `0` の級数を評価してもノルムは増えない |
| `hasEval_lubinTateEvalAtPoint` / `norm_lubinTateEvalAtPoint_lt_one` | 値も位相的冪零 |
| `eval_intertwine_mem_torsion` | `[π^n]_g(z) = 0` ならば `θ(z) ∈ Λ^f_n` |
| `coe_eval_at_eval_eq_eval_subst` | `u := η(z)` の**自身の座標系**での `θ(u)` は `z` の座標系での `(θ∘η)(z)` に等しい |
| `lubinTateClosure_le_of_intertwiner_section` | 上 2 つから `K_{π,f} ≤ K_{π,g}` |
| `lubinTateClosure_eq_of_same_uniformizer` | ★主定理 |

## 退化の自己検査

* **同じ `π` を落とすと偽**。`hf1 : coeff 1 f = π` と `hg1 : coeff 1 g = π`
  が同じ `π` であることは `LubinTateEndo`(`[1]_{g,f}` の存在)にも
  `powerSeries_uniqueness`(`π^{n+1} ≠ π` で割る)にも本質的に効いている。
  異なる素元 `π ≠ π'` では `[1]_{g,f}` が `𝒪_K` 係数では作れず、
  実際に別の塔になり得る(それが Λ6 本体)。
* **`f, g ∈ F_π` の条件を落とすと偽**。`coeff 1 = π` を落とすと
  `LubinTateEndo` の構成そのものが立たない。`map residue = X^q` を落とすと
  `Λ_n` の濃度・分離性(`iteratedLubinTateDistinguished` の Weierstrass 分解)
  が崩れ、`iteratedLubinTateTorsionPoints` の定義自体が意味を失う。
* **`coeff 1 θ = 1` を落とすと偽**。`a ≠ 1` の `[a]_{g,f}` は一般には
  全単射でない(`a` が単数でなければ捩れ点を潰す)。`a = 1` を使うので
  `θ ∘ η = η ∘ θ = X` が出る。

## 逸脱(記録)

* 原典 Milne Prop. 3.10 は `[u]_{g,f}` を**一般の単数 `u`** について作る。
  本ファイルは `u = 1` に限定する(Λ6a′ の担当範囲がそこだから)。
  一般の `u` は `K̂^ur` 係数・Dwork の補題が必要で、本ファイルの
  在庫だけでは届かない。★この限定により結論が弱まる箇所は無い
  (`K_{π,f} = K_{π,g}` は `u = 1` だけで出る)。
* 原典は「`K_{π,f}` と `K_{π,g}` の間の**同型**」ではなく、
  `K.closure` の中の**同じ中間体**であることを主張している。本ファイルも
  中間体の等式(`IntermediateField K.carrier K.closure` の `=`)で述べる。
* `lubinTateClosure_le_of_intertwiner_section` は `θ` について
  `constantCoeff θ = 0` も関数等式も**要求しない**——実際に効くのは
  `subst η θ = X`(`η` の切断であること)だけだと証明中に判明したため、
  仮定を落として一般化した。主定理の側では `θ := [1]_{g,f}` を渡すので
  実質同じことである。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 抽象核 —— 一般の可換環上の形式冪級数だけの主張 -/

/-- 定数項が `0` の冪級数の `n` 回自己合成も定数項が `0`。
`PowerSeries.constantCoeff_subst_eq_zero` による帰納法だけ。 -/
theorem constantCoeff_iteratedLubinTate_of_constantCoeff_zero {A : Type*} [CommRing A]
    {f : PowerSeries A} (hf0 : PowerSeries.constantCoeff f = 0) (n : ℕ) :
    PowerSeries.constantCoeff (iteratedLubinTate f n) = 0 := by
  induction n with
  | zero => exact PowerSeries.constantCoeff_X
  | succ m ih => exact PowerSeries.constantCoeff_subst_eq_zero ih f hf0

/-- 上を `PowerSeries.HasSubst` の形で言い直したもの。 -/
theorem hasSubst_iteratedLubinTate {A : Type*} [CommRing A]
    {f : PowerSeries A} (hf0 : PowerSeries.constantCoeff f = 0) (n : ℕ) :
    PowerSeries.HasSubst (iteratedLubinTate f n) := by
  show IsNilpotent (PowerSeries.constantCoeff (iteratedLubinTate f n))
  rw [constantCoeff_iteratedLubinTate_of_constantCoeff_zero hf0 n]
  exact IsNilpotent.zero

/-- ★★★★★★★★**抽象核 A**——`θ` が `f` と `g` を絡めるなら、その
`n` 回自己合成も絡める: `f ∘ θ = θ ∘ g` ならば `f^{(n)} ∘ θ = θ ∘ g^{(n)}`。

`PowerSeries.subst h α = α(h)` の規約で、`subst θ f = subst g θ` は
`f(θ(X)) = θ(g(X))` を意味する。証明は代入の結合律
(`PowerSeries.subst_comp_subst_apply`)を 3 回使う帰納法だけ——
`CommRing A` 以外に何も要らない(局所環も付値も分岐も出てこない)。 -/
theorem subst_iteratedLubinTate_of_intertwine {A : Type*} [CommRing A]
    {f g θ : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst θ f = PowerSeries.subst g θ) (n : ℕ) :
    PowerSeries.subst θ (iteratedLubinTate f n) =
      PowerSeries.subst (iteratedLubinTate g n) θ := by
  have hHSθ : PowerSeries.HasSubst θ := by
    show IsNilpotent (PowerSeries.constantCoeff θ); rw [hθ0]; exact IsNilpotent.zero
  have hHSg : PowerSeries.HasSubst g := by
    show IsNilpotent (PowerSeries.constantCoeff g); rw [hg0]; exact IsNilpotent.zero
  induction n with
  | zero =>
    show PowerSeries.subst θ PowerSeries.X = PowerSeries.subst PowerSeries.X θ
    rw [PowerSeries.subst_X hHSθ, PowerSeries.X_subst]
  | succ m ih =>
    show PowerSeries.subst θ (PowerSeries.subst (iteratedLubinTate f m) f) =
      PowerSeries.subst (PowerSeries.subst (iteratedLubinTate g m) g) θ
    rw [PowerSeries.subst_comp_subst_apply (hasSubst_iteratedLubinTate hf0 m) hHSθ f, ih,
      ← PowerSeries.subst_comp_subst_apply hHSθ (hasSubst_iteratedLubinTate hg0 m) f,
      hint,
      PowerSeries.subst_comp_subst_apply hHSg (hasSubst_iteratedLubinTate hg0 m) θ]

/-- ★★★★★★★★**抽象核 B**——互いに逆向きの絡み作用素で次数 1 の係数が
`1` のものは、合成が恒等になる: `f ∘ θ = θ ∘ g`・`g ∘ η = η ∘ f`・
`θ'(0) = η'(0) = 1` ならば `θ ∘ η = X`。

`ζ := θ ∘ η` が `f` との関数等式 `f ∘ ζ = ζ ∘ f` を満たし
(代入の結合律を 4 回)、次数 1 の係数が `1`
(`coeff_subst_eq_of_order_ge` の `n = 0` の場合)であることを示して、
一意性補題 `powerSeries_uniqueness` で `ζ = X` を出す。

分岐・Galois の語彙は 1 つも出てこないが、`powerSeries_uniqueness` が
`π^{n+1} - π` で割るために `IsLocalRing`・`IsDomain`・
`maximalIdeal = span {π}`・`π ≠ 0` を必要とする。 -/
theorem subst_eq_X_of_intertwine_pair {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A] {π : A}
    (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {f g θ η : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hg0 : PowerSeries.constantCoeff g = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ1 : PowerSeries.coeff 1 θ = 1)
    (hη0 : PowerSeries.constantCoeff η = 0) (hη1 : PowerSeries.coeff 1 η = 1)
    (hθint : PowerSeries.subst θ f = PowerSeries.subst g θ)
    (hηint : PowerSeries.subst η g = PowerSeries.subst f η) :
    PowerSeries.subst η θ = PowerSeries.X := by
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0]; exact IsNilpotent.zero
  have hHSg : PowerSeries.HasSubst g := by
    show IsNilpotent (PowerSeries.constantCoeff g); rw [hg0]; exact IsNilpotent.zero
  have hHSθ : PowerSeries.HasSubst θ := by
    show IsNilpotent (PowerSeries.constantCoeff θ); rw [hθ0]; exact IsNilpotent.zero
  have hHSη : PowerSeries.HasSubst η := by
    show IsNilpotent (PowerSeries.constantCoeff η); rw [hη0]; exact IsNilpotent.zero
  have hζ0 : PowerSeries.constantCoeff (PowerSeries.subst η θ) = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hη0 θ hθ0
  have hθord : ((0 + 1 : ℕ) : ℕ∞) ≤ θ.order := by
    simpa using PowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hθ0
  have hζ1 : PowerSeries.coeff 1 (PowerSeries.subst η θ) = 1 := by
    have := coeff_subst_eq_of_order_ge (δ := θ) (h := η) (π := (1 : A)) (n := 0) hη0 hη1 hθord hHSη
    simpa [hθ1] using this
  have hζfe : PowerSeries.subst f (PowerSeries.subst η θ) =
      PowerSeries.subst (PowerSeries.subst η θ) f := by
    rw [PowerSeries.subst_comp_subst_apply hHSη hHSf θ, ← hηint,
      ← PowerSeries.subst_comp_subst_apply hHSg hHSη θ, ← hθint,
      PowerSeries.subst_comp_subst_apply hHSθ hHSη f]
  exact powerSeries_uniqueness hπmax hπne0 hf0 hf1 hζ0 PowerSeries.constantCoeff_X
    (by rw [hζ1, PowerSeries.coeff_one_X]) hζfe
    (by rw [PowerSeries.subst_X hHSf, PowerSeries.X_subst])

/-! ## 1. 評価はノルムを増やさない(定数項 `0` の任意の冪級数について) -/

/-- **定数項 `0` の冪級数を評価してもノルムは増えない**: `‖φ(z)‖ ≤ ‖z‖`。
`Found/PGC/AdjoinIntegers.lean::norm_lubinTateActionAtTorsionPoint_le`
(`φ = [a]_f` に特化)と同じ `X` くくり出しの議論を、任意の `φ` へ
一般化しただけ。 -/
theorem norm_lubinTateEvalAtPoint_le {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (φ : PowerSeries (𝒪[K.carrier])) (hφ0 : PowerSeries.constantCoeff φ = 0) :
    ‖(↑(lubinTateEvalAtPoint K x z hz φ) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤
      ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr hφ0
  show ‖(↑(PowerSeries.aeval hz φ) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ _
  rw [hh, map_mul, aeval_X_eq_self]
  set y := PowerSeries.aeval hz h with hy_def
  show ‖(↑(z * y) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ _
  rw [Subring.coe_mul, norm_mul]
  calc ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ * ‖(↑y : _)‖
      ≤ ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ * 1 :=
        mul_le_mul_of_nonneg_left y.2 (norm_nonneg _)
    _ = _ := mul_one _

/-- 評価点のノルムが `1` 未満なら、値もまた位相的冪零——
`norm_lubinTateEvalAtPoint_le` と `hasEval_iff_coe` から。
これで値をさらに評価点として使える。 -/
theorem hasEval_lubinTateEvalAtPoint {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (hzlt : ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1)
    (φ : PowerSeries (𝒪[K.carrier])) (hφ0 : PowerSeries.constantCoeff φ = 0) :
    PowerSeries.HasEval (lubinTateEvalAtPoint K x z hz φ) := by
  rw [hasEval_iff_coe]
  exact tendsto_pow_atTop_nhds_zero_of_norm_lt_one
    (lt_of_le_of_lt (norm_lubinTateEvalAtPoint_le K x z hz φ hφ0) hzlt)

/-- 上と同じ状況で、値のノルムもまた `1` 未満。 -/
theorem norm_lubinTateEvalAtPoint_lt_one {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (hzlt : ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1)
    (φ : PowerSeries (𝒪[K.carrier])) (hφ0 : PowerSeries.constantCoeff φ = 0) :
    ‖(↑(lubinTateEvalAtPoint K x z hz φ) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1 :=
  lt_of_le_of_lt (norm_lubinTateEvalAtPoint_le K x z hz φ hφ0) hzlt

/-! ## 2. 絡み作用素は捩れ点を捩れ点へ送る -/

/-- ★★★★★★★★★★★★★★★★**絡み作用素は `g`-捩れ点を `f`-捩れ点へ送る**:
`[π^n]_g(z) = 0` かつ `f ∘ θ = θ ∘ g` ならば `θ(z) ∈ Λ^f_n`。

段取りは `Found/PGC/AdjoinIntegers.lean::lubinTateActionAtTorsionPoint_mem`
と全く同じ骨格:

1. 抽象核 A で `f^{(n)} ∘ θ = θ ∘ g^{(n)}`。
2. `aeval_subst_eq_aeval_aeval`(連鎖律)を 2 回使って
   `[π^n]_f(θ(z)) = θ([π^n]_g(z)) = θ(0) = 0`。
3. `eq_zero_of_pi_pow_action_eq_zero` で `D^f_n(θ(z)) = 0`。
4. `Polynomial.hom_eval₂` で `K.closure` のレベルへ押し出し、
   `iteratedLubinTateTorsionPoints` の定義に一致させる。

★`z` は `x` 自身である必要はない(`adjoinIntegers K x` の任意の
ノルム `< 1` の元でよい)。 -/
theorem eval_intertwine_mem_torsion {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (θ : PowerSeries (𝒪[K.carrier])) (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst θ f = PowerSeries.subst g θ)
    (n : ℕ) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (hzlt : ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1)
    (hzero : lubinTateEvalAtPoint K x z hz (LubinTateAction hq hπmax g hg0 hg1 hg (π ^ n)) = 0) :
    (↑(↑(lubinTateEvalAtPoint K x z hz θ) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hg0c : PowerSeries.constantCoeff g = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hg0
  have hHSθ : PowerSeries.HasSubst θ := by
    show IsNilpotent (PowerSeries.constantCoeff θ); rw [hθ0]; exact IsNilpotent.zero
  have hy : PowerSeries.HasEval (lubinTateEvalAtPoint K x z hz θ) :=
    hasEval_lubinTateEvalAtPoint K x z hz hzlt θ hθ0
  have hGn0 : PowerSeries.constantCoeff (iteratedLubinTate g n) = 0 :=
    constantCoeff_iteratedLubinTate_of_constantCoeff_zero hg0c n
  have hw : PowerSeries.HasEval (lubinTateEvalAtPoint K x z hz (iteratedLubinTate g n)) :=
    hasEval_lubinTateEvalAtPoint K x z hz hzlt _ hGn0
  have hweq : lubinTateEvalAtPoint K x z hz (iteratedLubinTate g n) = 0 := by
    rw [← LubinTateAction_pi_pow hq hπmax hπne0 g hg0 hg1 hg n]; exact hzero
  have hmain : lubinTateEvalAtPoint K x (lubinTateEvalAtPoint K x z hz θ) hy
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0 := by
    rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf n]
    show PowerSeries.aeval hy (iteratedLubinTate f n) = 0
    rw [← aeval_subst_eq_aeval_aeval (p := iteratedLubinTate f n) hHSθ hθ0 hz hy rfl,
      subst_iteratedLubinTate_of_intertwine hf0c hg0c hθ0 hint n,
      aeval_subst_eq_aeval_aeval (p := θ) (hasSubst_iteratedLubinTate hg0c n) hGn0 hz hw rfl]
    exact lubinTateEvalAtPoint_eq_zero_of_eq_zero K x hw hweq θ hθ0
  have hDn0 : Polynomial.aeval (lubinTateEvalAtPoint K x z hz θ)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 :=
    eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x _ hy hmain
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
  refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero, ?_⟩
  set G : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hG_def
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) G (lubinTateEvalAtPoint K x z hz θ)
  rw [← Polynomial.aeval_def] at key
  have hGcomp : G.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro u; rfl
  rw [hGcomp, hDn0, map_zero] at key
  show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).eval
      (G (lubinTateEvalAtPoint K x z hz θ)) = 0
  rw [Polynomial.eval_map]
  exact key.symm

/-! ## 3. 座標系をまたぐ評価(cross-point bridging) -/

/-- ★★★★★★★★★★★★**`u := η(z)` を自身の座標系で評価しても同じ値**:
`adjoinIntegers K u` の中で `θ` を `u` で評価したものは、
`adjoinIntegers K x` の中で `θ ∘ η` を `z` で評価したものに
(`K.closure` の値として)等しい。

`Found/PGC/LubinTateActionInclusion.lean::lubinTateEvalAtPoint_inclusion_comm`
(体の包含 `K⟮u⟯ ≤ K⟮x⟯` に沿った橋渡し)と
`aeval_subst_eq_aeval_aeval`(連鎖律)を繋ぐだけ。
★これが無いと「`θ(η(w)) = w`」が `w` の座標系でしか言えず、
`w ∈ K⟮η(w)⟯` が出ない。 -/
theorem coe_eval_at_eval_eq_eval_subst {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (hzlt : ‖(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1)
    (η : PowerSeries (𝒪[K.carrier])) (hη0 : PowerSeries.constantCoeff η = 0)
    (θ : PowerSeries (𝒪[K.carrier]))
    (u : K.closure)
    (hueq : u = (↑(↑(lubinTateEvalAtPoint K x z hz η) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure))
    (humem : u ∈ IntermediateField.adjoin K.carrier ({u} : Set K.closure))
    (huA : (⟨u, humem⟩ : IntermediateField.adjoin K.carrier ({u} : Set K.closure)) ∈
      adjoinIntegers K u)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({u} : Set K.closure))]
    (hu : PowerSeries.HasEval (⟨⟨u, humem⟩, huA⟩ : adjoinIntegers K u)) :
    (↑(↑(lubinTateEvalAtPoint K u ⟨⟨u, humem⟩, huA⟩ hu θ) :
        IntermediateField.adjoin K.carrier ({u} : Set K.closure)) : K.closure) =
      (↑(↑(lubinTateEvalAtPoint K x z hz (PowerSeries.subst η θ)) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hHSη : PowerSeries.HasSubst η := by
    show IsNilpotent (PowerSeries.constantCoeff η); rw [hη0]; exact IsNilpotent.zero
  have hle : IntermediateField.adjoin K.carrier ({u} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    apply IntermediateField.adjoin_simple_le_iff.mpr
    rw [hueq]; exact SetLike.coe_mem _
  have hincl : adjoinIntegersInclusionAlgHom K x u hle ⟨⟨u, humem⟩, huA⟩ =
      lubinTateEvalAtPoint K x z hz η := by
    apply Subtype.ext; apply Subtype.ext
    show (↑(↑((adjoinIntegersInclusion K x u hle) ⟨⟨u, humem⟩, huA⟩) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = _
    rw [coe_adjoinIntegersInclusion]
    exact hueq
  have hηeval : PowerSeries.HasEval (lubinTateEvalAtPoint K x z hz η) :=
    hasEval_lubinTateEvalAtPoint K x z hz hzlt η hη0
  have hu' : PowerSeries.HasEval (adjoinIntegersInclusionAlgHom K x u hle ⟨⟨u, humem⟩, huA⟩) := by
    rw [hincl]; exact hηeval
  have hcomm := lubinTateEvalAtPoint_inclusion_comm K x u hle ⟨⟨u, humem⟩, huA⟩ hu hu' θ
  rw [lubinTateEvalAtPoint_congr K x hu' hincl θ] at hcomm
  have hchain : lubinTateEvalAtPoint K x z hz (PowerSeries.subst η θ) =
      lubinTateEvalAtPoint K x (lubinTateEvalAtPoint K x z hz η) hηeval θ :=
    aeval_subst_eq_aeval_aeval (p := θ) hHSη hη0 hz hηeval rfl
  rw [hchain]
  exact hcomm.symm

/-! ## 4. `K_{π,f} ≤ K_{π,g}` と主定理 -/

/-- ★★★★★★★★★★★★★★★★**片側の包含**: `η` が `g` から `f` への絡み作用素
(`g ∘ η = η ∘ f`)で、`η` に切断 `θ`(`θ ∘ η = X`)があるなら
`K_{π,f} ≤ K_{π,g}`。

`w ∈ Λ^f_n` を取り、`u := η(w)`(`w` の座標系で計算)とすると
`eval_intertwine_mem_torsion` で `u ∈ Λ^g_n`、
`coe_eval_at_eval_eq_eval_subst` と `θ ∘ η = X` で
`θ(u) = w`(`u` の座標系で計算)、したがって `w ∈ K⟮u⟯ ≤ K_{π,g}`。

★`θ` については `subst η θ = X` 以外に何も要求しない(定数項が `0` で
あることも関数等式も、この補題の証明では使わない)。 -/
theorem lubinTateClosure_le_of_intertwiner_section {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (η : PowerSeries (𝒪[K.carrier])) (hη0 : PowerSeries.constantCoeff η = 0)
    (hηint : PowerSeries.subst η g = PowerSeries.subst f η)
    (θ : PowerSeries (𝒪[K.carrier])) (hcomp : PowerSeries.subst η θ = PowerSeries.X) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≤
      lubinTateClosure K hq hπmax hπne0 g hg0 hg1 hg := by
  rw [lubinTateClosure, IntermediateField.adjoin_le_iff]
  intro w hw
  simp only [lubinTateTorsionSet, Set.mem_iUnion, Finset.mem_coe] at hw
  obtain ⟨n, hwn⟩ := hw
  have hwmem : w ∈ IntermediateField.adjoin K.carrier ({w} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier w
  haveI : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({w} : Set K.closure)) :=
    finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n w hwn
  haveI := completeSpace_adjoinIntegers K w
  haveI := isLinearTopology_adjoinIntegers K w
  haveI := continuousSMul_adjoinIntegers K w
  set w' : adjoinIntegers K w := ⟨⟨w, hwmem⟩,
    mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n w hwn hwmem⟩ with hw'_def
  have hw' : PowerSeries.HasEval w' :=
    hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n w hwn hwmem
  have hwlt : ‖(↑w' : IntermediateField.adjoin K.carrier ({w} : Set K.closure))‖ < 1 := by
    show spectralNorm K.carrier K.closure w < 1
    exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n w hwn
  have hzero : lubinTateEvalAtPoint K w w' hw'
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0 :=
    pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n w hwn hwmem
  set u : K.closure := (↑(↑(lubinTateEvalAtPoint K w w' hw' η) :
    IntermediateField.adjoin K.carrier ({w} : Set K.closure)) : K.closure) with hu_def
  have hun : u ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 g hg0 hg1 hg n :=
    eval_intertwine_mem_torsion K hq hπmax hπne0 g hg0 hg1 hg f hf0 hf1 hf η hη0 hηint n w w' hw'
      hwlt hzero
  have humem : u ∈ IntermediateField.adjoin K.carrier ({u} : Set K.closure) :=
    IntermediateField.mem_adjoin_simple_self K.carrier u
  haveI : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({u} : Set K.closure)) :=
    finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 g hg0 hg1 hg n u hun
  have huA : (⟨u, humem⟩ : IntermediateField.adjoin K.carrier ({u} : Set K.closure)) ∈
      adjoinIntegers K u :=
    mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 g hg0 hg1 hg n u hun humem
  have hu : PowerSeries.HasEval (⟨⟨u, humem⟩, huA⟩ : adjoinIntegers K u) :=
    hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 g hg0 hg1 hg n u hun humem
  have hE := coe_eval_at_eval_eq_eval_subst K w w' hw' hwlt η hη0 θ u hu_def humem huA hu
  rw [hcomp, show lubinTateEvalAtPoint K w w' hw' PowerSeries.X = w' from aeval_X_eq_self hw'] at hE
  have hwK : w ∈ IntermediateField.adjoin K.carrier ({u} : Set K.closure) := by
    have h2 : (↑(↑(lubinTateEvalAtPoint K u ⟨⟨u, humem⟩, huA⟩ hu θ) :
        IntermediateField.adjoin K.carrier ({u} : Set K.closure)) : K.closure) = w := hE
    rw [← h2]
    exact SetLike.coe_mem _
  have huclos : u ∈ lubinTateClosure K hq hπmax hπne0 g hg0 hg1 hg := by
    apply IntermediateField.subset_adjoin
    simp only [lubinTateTorsionSet, Set.mem_iUnion, Finset.mem_coe]
    exact ⟨n, hun⟩
  exact (IntermediateField.adjoin_simple_le_iff.mpr huclos) hwK

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Λ6a′ ——`K_π` は `f` に依らない**: 同じ素元 `π` に属する 2 つの
Lubin-Tate 級数 `f, g ∈ F_π`(どちらも `coeff 0 = 0`・`coeff 1 = π`・
`map residue = X^q`)に対し `K_{π,f} = K_{π,g}`。

Milne CFT Prop. 3.10 の `u = 1` の場合。`θ := [1]_{g,f}`・`η := [1]_{f,g}`
(`LubinTateEndo`、どちらも `𝒪_K` 係数)を取り、抽象核 B で
`θ ∘ η = η ∘ θ = X` を出して、`lubinTateClosure_le_of_intertwiner_section`
を両向きに適用するだけ。★`K̂^ur`・`ℂ_K`・Dwork の補題はどれも通らない。

`π` そのものを取り替える場合(`π ≠ π'`)は本定理の射程外——それが
Λ6 本体(Lubin-Tate-Rosen)で、そちらは Dwork が要る。 -/
theorem lubinTateClosure_eq_of_same_uniformizer {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff)) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf =
      lubinTateClosure K hq hπmax hπne0 g hg0 hg1 hg := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hg0c : PowerSeries.constantCoeff g = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hg0
  set θ := LubinTateEndo hq hπmax g hg0c hg1 hg f hf0 hf1 hf 1 with hθ_def
  set η := LubinTateEndo hq hπmax f hf0c hf1 hf g hg0 hg1 hg 1 with hη_def
  have hθ0 : PowerSeries.constantCoeff θ = 0 :=
    constantCoeff_LubinTateEndo hq hπmax g hg0c hg1 hg f hf0 hf1 hf 1
  have hθ1 : PowerSeries.coeff 1 θ = 1 :=
    coeff_one_LubinTateEndo hq hπmax g hg0c hg1 hg f hf0 hf1 hf 1
  have hθint : PowerSeries.subst θ f = PowerSeries.subst g θ :=
    LubinTateEndo_functional_equation hq hπmax g hg0c hg1 hg f hf0 hf1 hf 1
  have hη0 : PowerSeries.constantCoeff η = 0 :=
    constantCoeff_LubinTateEndo hq hπmax f hf0c hf1 hf g hg0 hg1 hg 1
  have hη1 : PowerSeries.coeff 1 η = 1 :=
    coeff_one_LubinTateEndo hq hπmax f hf0c hf1 hf g hg0 hg1 hg 1
  have hηint : PowerSeries.subst η g = PowerSeries.subst f η :=
    LubinTateEndo_functional_equation hq hπmax f hf0c hf1 hf g hg0 hg1 hg 1
  have hcomp1 : PowerSeries.subst η θ = PowerSeries.X :=
    subst_eq_X_of_intertwine_pair hπmax hπne0 hf0c hf1 hg0c hθ0 hθ1 hη0 hη1 hθint hηint
  have hcomp2 : PowerSeries.subst θ η = PowerSeries.X :=
    subst_eq_X_of_intertwine_pair hπmax hπne0 hg0c hg1 hf0c hη0 hη1 hθ0 hθ1 hηint hθint
  exact le_antisymm
    (lubinTateClosure_le_of_intertwiner_section K hq hπmax hπne0 f hf0 hf1 hf g hg0 hg1 hg
      η hη0 hηint θ hcomp1)
    (lubinTateClosure_le_of_intertwiner_section K hq hπmax hπne0 g hg0 hg1 hg f hf0 hf1 hf
      θ hθ0 hθint η hcomp2)

end ABC3.Found.PGC
