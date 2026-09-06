import ABC3.Found.PGC.UnramifiedCompletionDVR
import ABC3.Found.PGC.ArithFrobeniusTopGen
import ABC3.Found.GaloisRep.AdicCompleteValued

/-!
# Dwork の補題(加法版)—— `σb - b = c` は `𝒪_{K̂^{ur}}` で常に解ける(`sorry` 無し)

経路 Λ の節点 **Λ6 の起点**。Λ6(`Art_π` の `π` 非依存性)の逐次近似は、
各段で「`σ` 差分方程式を解く」ことに帰着する。その加法版がこれである:

```
∃ σ ∈ Gal(K^ur/K), (σ は剰余体で z ↦ z^q) ∧ ∀ c ∈ 𝒪_{K̂^ur}, ∃ b ∈ 𝒪_{K̂^ur}, σb - b = c
```

★**量化の順が主張の中身**である。`∃ σ, ∀ c, ∃ b` であって `∀ c, ∃ σ, ∃ b` ではない
(後者は `c` ごとに Frobenius を選び直せるので空虚になる)。`σ` は `∃` の内側に
閉じ込めてあり、結論に自由なパラメータは出ない。

## 到達点

| | 主張 |
|---|---|
| 1 | `isAdicComplete_unramifiedCompletionInt`:**`𝒪_{K̂^{ur}}` は `𝔪` 進完備** |
| 2 | `exists_sub_frobenius_sub_eq_uniformizer_mul`:1 段の補正(剰余体で `t^q - t = c̄` を解く) |
| 3 | `dworkPartialSum_spec`:部分和の不変式 `c - (σb_N - b_N) = π^N c_N` |
| 4 | `exists_unramGalCompletionInt_sub_self_eq`:**`σb - b = c` の可解性**(`σ` は仮引数) |
| 5 | **`exists_arithFrobenius_dworkAdditive`**:主定理(`σ` を `∃` の内側に閉じ込めた形) |
| 6 | `exists_arithFrobenius_isCoherent_dworkAdditive`:同じ `σ` が位相的生成元でもある形 |

## ★見立てとの差(2026-09-06 の実測)

段取り係の見立ては「**解析ではなく純代数で書ける**——Cauchy 列・ε-δ を使わず
`IsPrecomplete.prec` + `IsHausdorff.haus` で済む」であった。★**当たっている**。
しかも**思ったより安い**:

* 段取り 1(「ここが最初の関門」と名指しされていた `IsAdicComplete` の instance)は
  **1 行で通った**。`isAdicComplete_valuationSubring` の側条件
  (`IsTopologicalRing` / `IsUniformAddGroup` on the subring)は**自動で付いた**——
  `unramifiedCompletionInt K` は `(Valued.v).valuationSubring` そのものなので
  部分環の位相・一様構造は mathlib の既製インスタンスが埋める。
  与える必要があったのは `isDiscreteValuationRing_unramifiedCompletionInt` だけで、
  それも `haveI` の型注釈(既定透明度)で `def` が展開されて通った。
* 逐次近似の不変式(段取り 5)は `linear_combination ih + π^N * h` の **1 行**で出た。
  帰納法の中で `σ` を分配するのは `map_add` / `map_mul` / `map_pow` と
  `unramGalCompletionInt_uniformizerCompletionInt`(`σπ = π`)の 4 つの `rw` だけ。
* ノルム・距離・収束は**一度も現れない**。本ファイルに `Metric` も `Tendsto` も無い。

★逆に高かったのは何も無い。全体で約 300 行、`lean_check` の往復は 10 回。

## 筋(7 段)

1. `IsAdicComplete (𝔪) 𝒪_{K̂^ur}`(`isAdicComplete_valuationSubring` に
   `CompleteSpace (K̂^ur)`(完備化なので自動)と DVR 性を与える)。
2. **1 段補題**:`residue (σe) = (residue e)^q` なので
   `residue (σe - e) = ē^q - ē`。剰余体は代数閉なので
   `exists_pow_card_sub_self_eq_completion` で `ē` が取れ、
   `IsLocalRing.residue_surjective` で `e` に持ち上がる。
3. 残差 `c - (σe - e)` は剰余 `0`、すなわち `𝔪 = (π)` の元なので `= π·c'`。
4. 再帰列 `c_0 := c`、`c_{n+1} := C(c_n)`、部分和 `b_N := Σ_{n<N} π^n E(c_n)`。
   ★`E` / `C` は `choose` で 1 度だけ取り出す(`Function.iterate` で回す)。
5. 不変式 `c - (σb_N - b_N) = π^N c_N`(`N` の帰納法 1 本)。要は **`σπ = π`**——
   `σ` は `K`-代数同型なので `algebraMap K.carrier _ π` を固定する。
6. `b_M - b_N ∈ (π^N)`(`N ≤ M`)から `IsPrecomplete.prec` で極限 `b` を取る。
7. `σb - b - c = σ(b - b_N) - (b - b_N) - (c - (σb_N - b_N))` の 3 項がすべて
   `(π^N)` に入るので、`IsHausdorff.haus` で `σb - b - c = 0`。

## 材料(すべて本プロジェクトの在庫)

| 要るもの | 出どころ |
|---|---|
| 算術 Frobenius(選択を 1 か所に閉じ込めた名前) | `ArithFrobeniusTopGen.arithFrobenius` / `arithFrobenius_residue` |
| `c^q - c = a` の可解性 | `UnramifiedResidueField.exists_pow_card_sub_self_eq_completion` |
| `𝔪 = (π)` / `π` で割る | `UnramifiedCompletionDVR.maximalIdeal_unramifiedCompletionInt_eq_span` |
| DVR 性 | `UnramifiedCompletionDVR.isDiscreteValuationRing_unramifiedCompletionInt` |
| `𝔪` 進完備性 | `GaloisRep.isAdicComplete_valuationSubring` |
| `σ` の `𝒪` への制限 | `UnramifiedCompletion.unramGalCompletionInt` |

★**`frobeniusCompletionInt` は使っていない**。あれは `coherentFrobenius`
(位相的生成元一般)から作られており、剰余体上の作用は `z ↦ z^{q^k}` でよい。
Dwork は**算術 Frobenius**でなければ回らない
(`UnramifiedResidueField.lean` の訂正欄)。

## ★設計上の注意(守ったこと)

* **既存の `Found/PGC/*.lean` を書き換えていない**(import のみ)。
* **`Skeleton` / `Interface` を触っていない**。`inertia` を経由せず、
  `Interface` の `SubgroupCorrespondence` / `ResidueCardinality` は
  主張にも証明にも現れない(Corollary 1.3 との循環を避けるため)。
* **`Abelianization` を使っていない**。
* **結論に自由なパラメータを出していない**——主定理の型は `K` にしか依存せず、
  `σ` も一意化元 `π` も `∃` / 証明の内側にある。
* **`sorry` 本体の `def` を作っていない**(D21)。

## 退化の自己検査

* **`σ` を任意の `unramGal K` に替えると偽**。`σ = 1` なら `σb - b = 0` なので
  `c ≠ 0` は解けない。★だから `σ` は `∃` の内側に無ければならない。
* **`σ` を「位相的生成元」に弱めるとこの証明は回らない**。剰余方程式が
  `ξ̄^{q^k} - ξ̄ = c̄`(`k ∈ Ẑ`)になり多項式方程式でなくなるので、段 2 の
  `exists_pow_card_sub_self_eq_completion` が当てられない。
  ★**主張自体が偽になるかは未確認**——「偽」とは書かない。
* `c ∈ 𝒪_{K̂^ur}` を落として `c ∈ K̂^{ur}` にしても主張は真だが(`π^{-n}` 倍で
  帰着する)、下流は `𝒪` 版しか使わないので示していない。
* `K^ur` を `K.closure` に替えると段 1 が**偽**——`𝒪_{ℂ_p}` の値群は稠密なので
  DVR でなく `𝔪` 進完備でもない。効いているのは不分岐性である。

## 逸脱(記録)

* **原典(Milne, CFT)の証明路を取っていない**。原典は `𝒪_{K̂^ur}` を
  `𝒪_{K^ur}` の `𝔪` 進完備化(逆極限)として作り、逆極限の各段で
  `A.7 / A.8` を回す。本ファイルはその路を**使わない**——Λ5 の逸脱
  (「ノルムの完備化と `𝔪` 進完備化の一致は未証明」)がまだ閉じていないためである。
  代わりに `IsAdicComplete` を**位相の側から**(完備な離散付値体の整数環は
  `𝔪` 進完備)取り、逐次近似はイデアルの言葉だけで回した。
  ★得られる主張は同じで、下流(Λ6 本体)への影響は無い。
* Λ5 の逸脱そのもの(ノルム完備化 ≅ `𝔪` 進完備化)は**依然として閉じていない**。
  本ファイルはそれを**必要としない**形に組み替えただけである。
* 乗法版(`σξ/ξ = u`、Dwork の補題の本来の形)は**入っていない**。本ファイルは
  加法版だけで、乗法版は別の節点(単数群の filtration が要る)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 0. `SModEq` とイデアルの往復

`IsPrecomplete` / `IsHausdorff` の形は `f m ≡ f n [SMOD I^m • ⊤]` である。
環を自分自身の加群と見ているので `I^m • ⊤ = I^m` で、差の所属に落とせる。 -/

/-- **`SModEq` を `I^n` への所属に落とす**——`I^n • (⊤ : Submodule R R) = I^n`。

★`IsPrecomplete.prec` / `IsHausdorff.haus` を使うたびに要る往復なので 1 度だけ書く。 -/
theorem sModEq_pow_iff {R : Type*} [CommRing R] (I : Ideal R) (n : ℕ) (x y : R) :
    x ≡ y [SMOD I ^ n • (⊤ : Submodule R R)] ↔ x - y ∈ I ^ n := by
  rw [SModEq.sub_mem]
  simp

/-! ## 1. `𝒪_{K̂^{ur}}` は `𝔪` 進完備

★段取りが「最初の関門」と名指しした所。実際には側条件が自動で付いた。 -/

/-- **★★★★★★★★★★★★(Λ6)`𝒪_{K̂^{ur}}` は `𝔪` 進完備**。

`GaloisRep.isAdicComplete_valuationSubring`(完備な離散付値体の整数環は
`𝔪` 進完備)に

* `CompleteSpace (K̂^{ur})`——`UniformSpace.Completion` なので自動、
* `IsDiscreteValuationRing 𝒪_{K̂^{ur}}`——`UnramifiedCompletionDVR` の到達点、
* `IsTopologicalRing` / `IsUniformAddGroup` on the subring——mathlib の
  部分環インスタンスが自動で埋める

を与えるだけ。★`unramifiedCompletionInt K` は `(Valued.v).valuationSubring` そのもの
なので、`haveI` の型注釈(既定透明度)で `def` が展開されて型が合う
(`instances` 透明度では合わないことに注意——`lean-idioms.md` の
「`↥(myDef)` に mathlib の `↥v.integer` 補題を `rw` できない」と同じ現象)。

退化の自己検査:`K^ur` を `K.closure` に替えると**偽**(`𝒪_{ℂ_p}` は DVR でない)。 -/
theorem isAdicComplete_unramifiedCompletionInt (K : PAdicLocalField p) :
    IsAdicComplete (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      ↥(unramifiedCompletionInt K) := by
  haveI : IsDiscreteValuationRing
      ((Valued.v : Valuation (unramifiedCompletion K) NNReal).valuationSubring) :=
    isDiscreteValuationRing_unramifiedCompletionInt K
  exact isAdicComplete_valuationSubring

/-- `𝔪^N = (π^N)`。逐次近似は `π` の冪の言葉で回すので、この形で持っておく。 -/
theorem maximalIdeal_pow_unramifiedCompletionInt (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (N : ℕ) :
    IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) ^ N
      = Ideal.span {uniformizerCompletionInt K π ^ N} := by
  rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax, Ideal.span_singleton_pow]

/-! ## 2. `σ` の基本性質 -/

/-- **算術 Frobenius の剰余体での作用(完備化側)**——`K^ur` 側の性質 `hσ` から出る。

★`UnramifiedResidueField.exists_arithFrobenius_completion` は `∃ σ` の形なので、
`arithFrobenius K`(`ArithFrobeniusTopGen`)という**特定の選択**には直接使えない。
同じ証明を `hσ` を仮引数に取る形で書き直す。 -/
theorem residue_unramGalCompletionInt (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (v : ↥(unramifiedCompletionInt K)) :
    IsLocalRing.residue _ (unramGalCompletionInt K σ v)
      = (IsLocalRing.residue _ v) ^ (Nat.card 𝓀[K.carrier]) := by
  obtain ⟨t, ht⟩ := (residueFieldEquiv K).surjective (IsLocalRing.residue _ v)
  calc IsLocalRing.residue _ (unramGalCompletionInt K σ v)
      = unramGalCompletionResidue K σ (IsLocalRing.residue _ v) :=
        (unramGalCompletionResidue_residue K σ v).symm
    _ = unramGalCompletionResidue K σ (residueFieldEquiv K t) := by rw [ht]
    _ = residueFieldEquiv K (unramGalResidue K σ t) :=
        (residueFieldEquiv_unramGalResidue K σ t).symm
    _ = residueFieldEquiv K (t ^ (Nat.card 𝓀[K.carrier])) := by
        obtain ⟨w, rfl⟩ := IsLocalRing.residue_surjective t
        rw [hσ w]
    _ = (residueFieldEquiv K t) ^ (Nat.card 𝓀[K.carrier]) := map_pow _ _ _
    _ = (IsLocalRing.residue _ v) ^ (Nat.card 𝓀[K.carrier]) := by rw [ht]

/-- **★証明の要:`σπ = π`**——`σ` は `K`-代数同型なので `algebraMap K.carrier _ π` を固定する。

これが無いと `§4` の不変式(`σ(π^n e) = π^n σe`)が崩れ、逐次近似が回らない。 -/
theorem unramGalCompletionInt_uniformizerCompletionInt (K : PAdicLocalField p) (σ : unramGal K)
    (π : 𝒪[K.carrier]) :
    unramGalCompletionInt K σ (uniformizerCompletionInt K π) = uniformizerCompletionInt K π := by
  apply Subtype.ext
  rw [unramGalCompletionInt_coe, uniformizerCompletionInt_coe]
  exact (unramGalCompletion K σ).commutes _

/-- `(π^N)` は `σ` で保たれる(`σπ = π` の系)。 -/
theorem unramGalCompletionInt_mem_span_pow (K : PAdicLocalField p) (σ : unramGal K)
    (π : 𝒪[K.carrier]) (N : ℕ) {x : ↥(unramifiedCompletionInt K)}
    (hx : x ∈ Ideal.span {uniformizerCompletionInt K π ^ N}) :
    unramGalCompletionInt K σ x ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
  rw [Ideal.mem_span_singleton'] at hx ⊢
  obtain ⟨a, ha⟩ := hx
  exact ⟨unramGalCompletionInt K σ a, by
    rw [← ha, map_mul, map_pow, unramGalCompletionInt_uniformizerCompletionInt]⟩

/-! ## 3. 1 段の補正 -/

/-- **★★★★★★★★★★(Λ6)1 段の補正**——任意の `c ∈ 𝒪_{K̂^{ur}}` に対し
`c - (σe - e) = π·c'` なる `e`, `c'` が取れる。

剰余体で `t^q - t = c̄` を解き(`exists_pow_card_sub_self_eq_completion`、
`𝓀_{K̂^{ur}}` が代数閉であることによる)、`residue_surjective` で `e` に持ち上げる。
`residue (σe) = ē^q` なので残差 `c - (σe - e)` は剰余 `0`、すなわち `𝔪 = (π)` の元。

退化の自己検査。

* `hσ`(`σ` が剰余体で `q` 乗)を落とすと**偽**——`σ = 1` なら
  `c - 0 = c` は一般に `𝔪` に入らない。
* `hπmax` を落とすと `π` が自由なパラメータになり、`π` が単元のとき主張は空虚
  (逆に `π = π₀²` では `c'` が取れない)。 -/
theorem exists_sub_frobenius_sub_eq_uniformizer_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (c : ↥(unramifiedCompletionInt K)) :
    ∃ e : ↥(unramifiedCompletionInt K), ∃ c' : ↥(unramifiedCompletionInt K),
      c - (unramGalCompletionInt K σ e - e) = uniformizerCompletionInt K π * c' := by
  obtain ⟨t, ht⟩ := exists_pow_card_sub_self_eq_completion K (IsLocalRing.residue _ c)
  obtain ⟨e, rfl⟩ := IsLocalRing.residue_surjective t
  refine ⟨e, ?_⟩
  have hmem : c - (unramGalCompletionInt K σ e - e)
      ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
    rw [← IsLocalRing.residue_eq_zero_iff, map_sub, map_sub,
      residue_unramGalCompletionInt K hσ e, ht]
    ring
  rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax,
    Ideal.mem_span_singleton'] at hmem
  obtain ⟨c', hc'⟩ := hmem
  exact ⟨c', by rw [← hc', mul_comm]⟩

/-! ## 4. 逐次近似の部分和とその不変式

`E` / `C` は `§3` を `choose` して得る 2 つの写像で、
`c_n := C^[n] c`、`e_n := E (C^[n] c)`、`b_N := Σ_{n<N} π^n e_n` である。 -/

/-- 逐次近似の部分和 `b_N = Σ_{n<N} π^n E(C^[n] c)`。 -/
noncomputable def dworkPartialSum (K : PAdicLocalField p) (π : 𝒪[K.carrier])
    (E C : ↥(unramifiedCompletionInt K) → ↥(unramifiedCompletionInt K))
    (c : ↥(unramifiedCompletionInt K)) (N : ℕ) : ↥(unramifiedCompletionInt K) :=
  ∑ n ∈ Finset.range N, uniformizerCompletionInt K π ^ n * E (C^[n] c)

theorem dworkPartialSum_zero (K : PAdicLocalField p) (π : 𝒪[K.carrier])
    (E C : ↥(unramifiedCompletionInt K) → ↥(unramifiedCompletionInt K))
    (c : ↥(unramifiedCompletionInt K)) : dworkPartialSum K π E C c 0 = 0 := by
  rw [dworkPartialSum, Finset.range_zero, Finset.sum_empty]

theorem dworkPartialSum_succ (K : PAdicLocalField p) (π : 𝒪[K.carrier])
    (E C : ↥(unramifiedCompletionInt K) → ↥(unramifiedCompletionInt K))
    (c : ↥(unramifiedCompletionInt K)) (N : ℕ) :
    dworkPartialSum K π E C c (N + 1)
      = dworkPartialSum K π E C c N + uniformizerCompletionInt K π ^ N * E (C^[N] c) := by
  rw [dworkPartialSum, dworkPartialSum, Finset.sum_range_succ]

/-- **★★★★★★★★★★★★(Λ6)逐次近似の不変式**:`c - (σb_N - b_N) = π^N c_N`。

`N` の帰納法 1 本。`σ` を和と積に分配するのに `map_add` / `map_mul` / `map_pow` と
**`σπ = π`**(`unramGalCompletionInt_uniformizerCompletionInt`)を使い、
あとは `linear_combination ih + π^N * h` で閉じる。

退化の自己検査:`σπ = π` を落とすと `σ(π^n e) = π^n σe` が言えず、帰納段が崩れる。 -/
theorem dworkPartialSum_spec (K : PAdicLocalField p) (σ : unramGal K) {π : 𝒪[K.carrier]}
    {E C : ↥(unramifiedCompletionInt K) → ↥(unramifiedCompletionInt K)}
    (hEC : ∀ x, x - (unramGalCompletionInt K σ (E x) - E x)
      = uniformizerCompletionInt K π * C x)
    (c : ↥(unramifiedCompletionInt K)) (N : ℕ) :
    c - (unramGalCompletionInt K σ (dworkPartialSum K π E C c N)
        - dworkPartialSum K π E C c N)
      = uniformizerCompletionInt K π ^ N * C^[N] c := by
  induction N with
  | zero =>
    rw [dworkPartialSum_zero, map_zero, pow_zero, one_mul, Function.iterate_zero_apply]
    ring
  | succ N ih =>
    have h := hEC (C^[N] c)
    rw [dworkPartialSum_succ, map_add, map_mul, map_pow,
      unramGalCompletionInt_uniformizerCompletionInt, Function.iterate_succ_apply', pow_succ]
    linear_combination ih + uniformizerCompletionInt K π ^ N * h

/-- **部分和は `(π^N)` を法として安定**:`b_{N+k} - b_N ∈ (π^N)`。

★`IsPrecomplete.prec` が要求する Cauchy 条件はこれ(と差の符号の入れ替え)である。 -/
theorem dworkPartialSum_sub_mem (K : PAdicLocalField p) (π : 𝒪[K.carrier])
    (E C : ↥(unramifiedCompletionInt K) → ↥(unramifiedCompletionInt K))
    (c : ↥(unramifiedCompletionInt K)) (N k : ℕ) :
    dworkPartialSum K π E C c (N + k) - dworkPartialSum K π E C c N
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
  induction k with
  | zero => simp
  | succ k ih =>
    have hstep : dworkPartialSum K π E C c (N + (k + 1))
        = dworkPartialSum K π E C c (N + k)
          + uniformizerCompletionInt K π ^ (N + k) * E (C^[N + k] c) := by
      rw [← Nat.add_assoc, dworkPartialSum_succ]
    have hmem : uniformizerCompletionInt K π ^ (N + k) * E (C^[N + k] c)
        ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
      rw [Ideal.mem_span_singleton']
      exact ⟨uniformizerCompletionInt K π ^ k * E (C^[N + k] c), by rw [pow_add]; ring⟩
    have heq : dworkPartialSum K π E C c (N + (k + 1)) - dworkPartialSum K π E C c N
        = (dworkPartialSum K π E C c (N + k) - dworkPartialSum K π E C c N)
          + uniformizerCompletionInt K π ^ (N + k) * E (C^[N + k] c) := by
      rw [hstep]; ring
    rw [heq]
    exact Submodule.add_mem _ ih hmem

/-! ## 5. 可解性 -/

/-- **★★★★★★★★★★★★★★★★★★(Λ6)`σb - b = c` は `𝒪_{K̂^{ur}}` で解ける**
(`σ` を仮引数に取った形)。

`§3` を `choose` して `E` / `C` を取り、部分和 `b_N`(`§4`)を作る。
`§4` の不変式と `dworkPartialSum_sub_mem` から `IsPrecomplete.prec` で極限 `b` を取り、

```
σb - b - c = σ(b - b_N) - (b - b_N) - (c - (σb_N - b_N))
```

の 3 項がすべて `(π^N)` に入ることを見て `IsHausdorff.haus` で `0` にする。

★ノルム・距離・収束は使わない。`𝔪` 進完備性(`§1`)だけで閉じる。

退化の自己検査:`hσ` を落とすと**偽**(`σ = 1` は `hσ` を満たさず、
`σb - b = 0` なので `c ≠ 0` は解けない)。 -/
theorem exists_unramGalCompletionInt_sub_self_eq (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (c : ↥(unramifiedCompletionInt K)) :
    ∃ b : ↥(unramifiedCompletionInt K), unramGalCompletionInt K σ b - b = c := by
  classical
  haveI hAC : IsAdicComplete (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      ↥(unramifiedCompletionInt K) := isAdicComplete_unramifiedCompletionInt K
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  choose E C hEC using exists_sub_frobenius_sub_eq_uniformizer_mul K hσ hπ0 hπmax
  have hpow := maximalIdeal_pow_unramifiedCompletionInt K hπ0 hπmax
  have hcauchy : ∀ {m n : ℕ}, m ≤ n →
      dworkPartialSum K π E C c m ≡ dworkPartialSum K π E C c n
        [SMOD (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)) ^ m
          • (⊤ : Submodule ↥(unramifiedCompletionInt K) ↥(unramifiedCompletionInt K))] := by
    intro m n hmn
    rw [sModEq_pow_iff, hpow]
    have h := dworkPartialSum_sub_mem K π E C c m (n - m)
    rw [Nat.add_sub_cancel' hmn] at h
    simpa using Submodule.neg_mem _ h
  obtain ⟨b, hb⟩ := IsPrecomplete.prec inferInstance hcauchy
  refine ⟨b, sub_eq_zero.mp ?_⟩
  refine IsHausdorff.haus (I := IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
    inferInstance _ (fun N => ?_)
  rw [sModEq_pow_iff, sub_zero, hpow]
  have h1 : b - dworkPartialSum K π E C c N
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
    have h := (sModEq_pow_iff _ N _ _).mp (hb N)
    rw [hpow] at h
    simpa using Submodule.neg_mem _ h
  have h2 : c - (unramGalCompletionInt K σ (dworkPartialSum K π E C c N)
      - dworkPartialSum K π E C c N) ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
    rw [dworkPartialSum_spec K σ hEC c N, Ideal.mem_span_singleton']
    exact ⟨C^[N] c, mul_comm _ _⟩
  have h3 := unramGalCompletionInt_mem_span_pow K σ π N h1
  have hdecomp : unramGalCompletionInt K σ b - b - c
      = (unramGalCompletionInt K σ (b - dworkPartialSum K π E C c N))
        - (b - dworkPartialSum K π E C c N)
        - (c - (unramGalCompletionInt K σ (dworkPartialSum K π E C c N)
            - dworkPartialSum K π E C c N)) := by
    rw [map_sub]; ring
  rw [hdecomp]
  exact Submodule.sub_mem _ (Submodule.sub_mem _ h3 h1) h2

/-- `arithFrobenius K`(選択を 1 か所に閉じ込めた算術 Frobenius)版。下流はこれを使う。 -/
theorem exists_arithFrobenius_sub_self_eq (K : PAdicLocalField p)
    (c : ↥(unramifiedCompletionInt K)) :
    ∃ b : ↥(unramifiedCompletionInt K),
      unramGalCompletionInt K (arithFrobenius K) b - b = c :=
  exists_unramGalCompletionInt_sub_self_eq K (arithFrobenius_residue K) c

/-- `b ↦ σb - b` は `𝒪_{K̂^{ur}}` 上全射(算術 Frobenius について)。 -/
theorem surjective_unramGalCompletionInt_sub_self (K : PAdicLocalField p) :
    Function.Surjective
      (fun b : ↥(unramifiedCompletionInt K) =>
        unramGalCompletionInt K (arithFrobenius K) b - b) :=
  fun c => exists_arithFrobenius_sub_self_eq K c

/-! ## 6. 主定理 -/

/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ6 の起点)Dwork の補題(加法版)**——
`Gal(K^ur/K)` には、

* `𝓀_{K^{ur}}` 上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
* **その 1 つの `σ` について**、任意の `c ∈ 𝒪_{K̂^{ur}}` に対し
  `σb - b = c` なる `b ∈ 𝒪_{K̂^{ur}}` が存在する

元 `σ` がある。

★**量化の順が主張の中身**である。`∃ σ, ∀ c, ∃ b` を主張しており、
`∀ c, ∃ σ, ∃ b` ではない(後者は `c` ごとに Frobenius を選び直せるので空虚)。

**証明**:`arithFrobenius K`(`ArithFrobeniusTopGen`)に
`exists_unramGalCompletionInt_sub_self_eq` を当てる。

退化の自己検査。

* **`σ` を任意の `unramGal K` に替えると偽**——`σ = 1` なら `σb - b = 0`。
* `σ` を「位相的生成元」に弱めると**この証明は回らない**(剰余方程式が
  `ξ̄^{q^k} - ξ̄ = c̄` になり多項式方程式でなくなる)。
  ★主張自体が偽になるかは**未確認**。
* `c ∈ 𝒪_{K̂^{ur}}` を落としても主張は真だが、下流は `𝒪` 版しか使わない。 -/
theorem exists_arithFrobenius_dworkAdditive (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      ∀ c : ↥(unramifiedCompletionInt K), ∃ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b - b = c :=
  ⟨arithFrobenius K, arithFrobenius_residue K, fun c => exists_arithFrobenius_sub_self_eq K c⟩

/-- **★★★★★★★★★★★★★★★★★★★★同じ `σ` が位相的生成元でもある形**。

Λ5(`Ẑ` との同定)と Λ6(Dwork)が**同じ 1 つの元**を指せるように、
`ArithFrobeniusTopGen.exists_arithFrobenius_isCoherent` の 2 性質に
可解性を並べた形も出しておく。★3 つを別々の `∃` に分けたら無意味になる。 -/
theorem exists_arithFrobenius_isCoherent_dworkAdditive (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
      ∀ c : ↥(unramifiedCompletionInt K), ∃ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b - b = c :=
  ⟨arithFrobenius K, arithFrobenius_residue K, arithFrobenius_mem_unramLevelGeneratorSet K,
    fun c => exists_arithFrobenius_sub_self_eq K c⟩

end ABC3.Found.PGC
