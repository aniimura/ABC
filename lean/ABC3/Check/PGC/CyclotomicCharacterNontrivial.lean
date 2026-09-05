import ABC3.Skeleton.PGC.Section1
import ABC3.Found.PGC.QpRootsOfUnity
import Mathlib.FieldTheory.Galois.Infinite
import Mathlib.RingTheory.RootsOfUnity.AlgebraicallyClosed
import Mathlib.Data.Nat.Choose.Dvd

/-!
# `χ_{ℚ_p}` は自明指標ではない——Proposition 1.1 の非退化検査

`Skeleton/PGC/Section1.lean` の `cyclotomicCharacterObject` は

```
obj := fun K g => cyclotomicCharacter K.closure p g.toRingEquiv
```

と定義されている。ところが mathlib の `cyclotomicCharacter` の docstring は

> Note: This is the trivial character when `L` does not contain all `pⁱ`-th primitive roots.

と明記している——**十分な 1 の冪根が無ければ定数 1 として定義される**。
もし `obj K` が定数関数 `1` なら、Proposition 1.1(`RecoverableFromAbsGal`、
「α による移送で対応が保たれる」)は**両辺とも定数 1** になり、内容ゼロで真になる。

この木では「落とした条件は主張を偽にするか自明になるかのどちらかになる」例が
6 つ見つかっている(`Check/PGC/*Degenerate.lean`)。Proposition 1.1 に
残っていた可能性は**後者**——`Check/PGC/RefutationAttempts.lean` は
「共役による経路は閉じている」ことしか確かめておらず、
「そもそも χ が定数 1 ではないか」は検査していなかった。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。

回復されると言われているのは、まさにこの χ である。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

## 何を示したか

`cyclotomicCharacter_selfField_ne_one`:
`∃ g : (selfField p).absGal, cyclotomicCharacter (selfField p).closure p g.toRingEquiv ≠ 1`
——`sorry` 無し、すべての素数 `p`(`p = 2` を含む)について。

したがって `cyclotomicCharacterObject.obj (selfField p)` は定数関数 1 ではなく
(`cyclotomicCharacterObject_selfField_ne_const_one`)、Proposition 1.1 は
少なくとも `K = ℚ_p` において**空でない内容**を主張している。

## 証明の筋

1. `AlgebraicClosure.hasEnoughRootsOfUnity_pow` により `K̄ = closure ℚ_p` は
   すべての `p^i` 乗根を持つので、`cyclotomicCharacter.spec` が使える:
   `t ^ p ^ n = 1` なら `g t = t ^ ((χ g).val.toZModPow n).val`。
2. χ ≡ 1 と仮定すると指数は `(1 : ZMod (p^n)).val = 1` なので、
   **すべての p 冪乗根が Γ_{ℚ_p} の全元で固定される**。
3. `K̄/ℚ_p` は Galois(標数 0)なので `fixedField ⊤ = ⊥`
   (`InfiniteGalois.fixedField_fixingSubgroup` を `⊥` に適用)。
   よって原始 `p²` 乗根 ζ が `ℚ_p` の中に入る。
4. これは偽である。`n = 2` に取ったのは `p = 2` を込みで一様に扱うため:
   * `p ≥ 3` なら `y := ζ^p` は `y^p = 1`・`y ≠ 1` を満たす原始 `p` 乗根。
   * `p = 2` なら ζ 自身が原始 4 乗根で `ζ² ≠ 1`。
5. `ℚ_p` に原始 `p` 乗根が無いこと(`p ≥ 3`)は
   **ノルムを使わない整除の議論**で出る(`eq_zero_of_add_one_pow_eq_one`):
   `b := ζ - 1` は `p ∣ b`(剰余体 `ZMod p` で `ZMod.pow_card` を使う)を満たし、
   二項展開 `(1+b)^p = 1` の `k ≥ 2` の項はすべて `p·b²` で割れるので
   `p·b² ∣ p·b`、すなわち `b² ∣ b`。`b ≠ 0` なら `b` は単元となり、
   `p ∣ b` から `p` が `ℤ_p` の単元になって矛盾。
   ★`p = 2` ではこの議論は**破れる**(ζ = −1 が実際に `ℚ_2` にある)。
6. `p = 2` は `ζ² = −1` を `ZMod 4` に落として `decide` で潰す。

## 逸脱

無し。原典の主張についてではなく、我々の形式化が退化していないことの検査である。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC Finset

/-! ## `ℚ_p` に原始 `p` 乗根は無い(`p ≥ 3`) -/

/-- `∑_{k<n+1} f k` から `k = 0, 1` の 2 項を剥がす。
`n` の側に依存型が乗らないよう `n = m + 2` に分解してから使う。 -/
private theorem sum_range_split_two {M : Type*} [AddCommMonoid M] {n : ℕ} (hn : 2 ≤ n)
    (f : ℕ → M) :
    ∑ k ∈ range (n + 1), f k = (∑ k ∈ range (n - 1), f (k + 2)) + f 1 + f 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  rw [Finset.sum_range_succ' f (m + 2)]
  congr 1
  rw [Finset.sum_range_succ' (fun k => f (k + 1)) (m + 1)]
  simp

/-- **`p ≥ 3` なら `(1 + b)^p = 1` かつ `p ∣ b` を満たす `b : ℤ_[p]` は `0` だけ。**

二項展開の `k ≥ 2` の項はすべて `p · b²` で割れる:
`2 ≤ k < p` では `p ∣ C(p,k)` かつ `b² ∣ b^k`、`k = p` では
`b^p = b² · b^{p-2}` で `p ∣ b ∣ b^{p-2}`(ここで `p ≥ 3` を使う)。
`k = 1` の項が `p · b` なので `p · b² ∣ p · b`、すなわち `b² ∣ b`。
`b ≠ 0` を消去すると `b ∣ 1` となり、`p ∣ b` から `p` が単元——矛盾。

★`p = 2` では `k = p = 2` の項が `b²` で、`p ∣ b^{p-2} = b^0 = 1` が偽になる。
実際 `ζ = -1`(`b = -2`)が反例なので、この仮定は落とせない。 -/
theorem eq_zero_of_add_one_pow_eq_one (p : ℕ) [hp : Fact p.Prime] (hp3 : 3 ≤ p) (b : ℤ_[p])
    (hdvd : (p : ℤ_[p]) ∣ b) (h : (b + 1) ^ p = 1) : b = 0 := by
  by_contra hb0
  set f : ℕ → ℤ_[p] := fun k => b ^ k * (p.choose k : ℤ_[p]) with hf
  have hexp : (b + 1) ^ p = ∑ k ∈ range (p + 1), f k := by
    rw [add_pow]; simp [hf]
  have hf0 : f 0 = 1 := by simp [hf]
  have hf1 : f 1 = b * (p : ℤ_[p]) := by simp [hf]
  have key : (∑ k ∈ range (p - 1), f (k + 2)) + b * (p : ℤ_[p]) = 0 := by
    have h2 := hexp.symm.trans h
    rw [sum_range_split_two (by omega) f, hf0, hf1] at h2
    have h3 : (∑ k ∈ range (p - 1), f (k + 2)) + b * (p : ℤ_[p]) + 1 = 0 + 1 := by
      rw [zero_add]; exact h2
    exact add_right_cancel h3
  have hdvdS : ((p : ℤ_[p]) * b ^ 2) ∣ (∑ k ∈ range (p - 1), f (k + 2)) := by
    refine Finset.dvd_sum ?_
    intro k hk
    rw [Finset.mem_range] at hk
    simp only [hf]
    rcases Nat.lt_or_ge (k + 2) p with hlt | hge
    · have hc : (p : ℤ_[p]) ∣ (p.choose (k + 2) : ℤ_[p]) :=
        Nat.cast_dvd_cast (hp.out.dvd_choose_self (by omega) hlt)
      calc (p : ℤ_[p]) * b ^ 2 ∣ (p.choose (k + 2) : ℤ_[p]) * b ^ (k + 2) :=
            mul_dvd_mul hc (pow_dvd_pow b (by omega))
        _ = b ^ (k + 2) * (p.choose (k + 2) : ℤ_[p]) := mul_comm _ _
    · have heq : k + 2 = p := by omega
      rw [heq, Nat.choose_self, Nat.cast_one, mul_one]
      have hsp : b ^ p = b ^ 2 * b ^ (p - 2) := by rw [← pow_add]; congr 1; omega
      rw [hsp, mul_comm ((p : ℤ_[p])) (b ^ 2)]
      exact mul_dvd_mul_left _ (dvd_trans hdvd (dvd_pow_self b (by omega)))
  have hpb : ((p : ℤ_[p]) * b ^ 2) ∣ ((p : ℤ_[p]) * b) := by
    have hng : (p : ℤ_[p]) * b = -(∑ k ∈ range (p - 1), f (k + 2)) := by linear_combination key
    rw [hng]
    exact dvd_neg.mpr hdvdS
  have hpne : ((p : ℤ_[p]) ≠ 0) := by exact_mod_cast (NeZero.ne (p : ℤ_[p]))
  have hb2b : b ^ 2 ∣ b := (mul_dvd_mul_iff_left hpne).mp hpb
  have hbb : b * b ∣ b * 1 := by rwa [mul_one, ← pow_two]
  have hb1 : b ∣ 1 := (mul_dvd_mul_iff_left hb0).mp hbb
  have hpu : IsUnit ((p : ℤ_[p])) := isUnit_of_dvd_one (hdvd.trans hb1)
  rw [PadicInt.isUnit_iff, PadicInt.norm_p, inv_eq_one] at hpu
  have h1 : (1 : ℝ) < p := by exact_mod_cast hp.out.one_lt
  linarith

/-- **`p ≥ 3` のとき `ℚ_p` の中の `p` 乗根は `1` だけ。**

`ℚ_p` に原始 `p` 乗根が入らないという、円分指標の非自明性の核。
`Found/PGC/QpRootsOfUnity.lean::eq_one_of_pow_eq_one_qp` は `ℓ ≠ p` の場合
(`p` と素な捩れの計数)なので、`ℓ = p` の本件は別立てになる。 -/
theorem eq_one_of_pow_p_eq_one_qp (p : ℕ) [hp : Fact p.Prime] (hp3 : 3 ≤ p) {ζ : ℚ_[p]}
    (h : ζ ^ p = 1) : ζ = 1 := by
  have hnorm : ‖ζ‖ = 1 :=
    eq_one_of_pow_eq_one_real (x := ‖ζ‖) (n := p) hp.out.ne_zero (norm_nonneg _)
      (by rw [← norm_pow, h, norm_one])
  set z : ℤ_[p] := ⟨ζ, le_of_eq hnorm⟩ with hz
  have hzp : z ^ p = 1 := by
    apply Subtype.ext
    show ζ ^ p = (1 : ℚ_[p])
    exact h
  have hdvd : (p : ℤ_[p]) ∣ (z - 1) := by
    have h1 : PadicInt.toZMod z = 1 := by
      have hc := congrArg PadicInt.toZMod hzp
      rw [map_pow, map_one, ZMod.pow_card] at hc
      exact hc
    have h2 : PadicInt.toZMod (z - 1) = 0 := by rw [map_sub, h1, map_one, sub_self]
    have h3 : (z - 1) ∈ RingHom.ker (PadicInt.toZMod (p := p)) := h2
    rw [PadicInt.ker_toZMod, PadicInt.maximalIdeal_eq_span_p, Ideal.mem_span_singleton] at h3
    exact h3
  have hb : z - 1 = 0 :=
    eq_zero_of_add_one_pow_eq_one p hp3 (z - 1) hdvd (by rw [sub_add_cancel]; exact hzp)
  have hz1 : z = 1 := by linear_combination hb
  have hfin := congrArg (fun w : ℤ_[p] => (w : ℚ_[p])) hz1
  simpa [hz] using hfin

/-- **`ℚ_2` の中の 4 乗根は 2 乗根。**——`ζ² = -1` を `ZMod 4` に落とすと
平方が `{0, 1}` にしかならないので潰れる(`-1 = 3`)。
すなわち `ℚ_2` に原始 4 乗根 `i` は無い。 -/
theorem sq_eq_one_of_pow_four_eq_one_q2 {ζ : ℚ_[2]} (h : ζ ^ 4 = 1) : ζ ^ 2 = 1 := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  have hnorm : ‖ζ‖ = 1 :=
    eq_one_of_pow_eq_one_real (x := ‖ζ‖) (n := 4) (by norm_num) (norm_nonneg _)
      (by rw [← norm_pow, h, norm_one])
  set z : ℤ_[2] := ⟨ζ, le_of_eq hnorm⟩ with hz
  have hz4 : z ^ 4 = 1 := by
    apply Subtype.ext
    show ζ ^ 4 = (1 : ℚ_[2])
    exact h
  have hfac : (z ^ 2 - 1) * (z ^ 2 + 1) = 0 := by linear_combination hz4
  rcases mul_eq_zero.mp hfac with h1 | h2
  · have hsq : z ^ 2 = 1 := by linear_combination h1
    have hfin := congrArg (fun w : ℤ_[2] => (w : ℚ_[2])) hsq
    simpa [hz] using hfin
  · exfalso
    have h3 : z ^ 2 = -1 := by linear_combination h2
    have h4 := congrArg (PadicInt.toZModPow (p := 2) 2) h3
    rw [map_pow, map_neg, map_one] at h4
    exact (by decide : ∀ x : ZMod (2 ^ 2), x ^ 2 ≠ -1) _ h4

/-! ## 本題: χ は定数 1 ではない -/

/-- **★★★★★★円分指標 `χ : Γ_{ℚ_p} → ℤ_[p]ˣ` は自明指標ではない。**

mathlib の `cyclotomicCharacter` は「`L` が十分な `pⁱ` 乗根を持たないときは
定数 1」という規約で定義されている。`Skeleton/PGC/Section1.lean` の
`cyclotomicCharacterObject` がその規約に落ちていないこと——すなわち
Proposition 1.1 が内容ゼロで真になっていないこと——を、`K = ℚ_p` で確かめる。

証明はモジュール docstring の 1〜6 のとおり。 -/
theorem cyclotomicCharacter_selfField_ne_one (p : ℕ) [hp : Fact p.Prime] :
    ∃ g : (selfField p).absGal,
      cyclotomicCharacter (selfField p).closure p g.toRingEquiv ≠ 1 := by
  by_contra hcon
  simp only [not_exists, not_not] at hcon
  -- 1. `K̄` は全ての `p^i` 乗根を持つ
  haveI hEnough : ∀ i : ℕ, HasEnoughRootsOfUnity (selfField p).closure (p ^ i) :=
    fun i => AlgebraicClosure.hasEnoughRootsOfUnity_pow ℚ_[p] p i
  obtain ⟨ζ, hζ⟩ := HasEnoughRootsOfUnity.exists_primitiveRoot (selfField p).closure (p ^ 2)
  have hζpow : ζ ^ p ^ 2 = 1 := hζ.pow_eq_one
  haveI : Fact (1 < p ^ 2) := ⟨by have := hp.out.two_le; nlinarith⟩
  -- 2. χ ≡ 1 なら指数は 1 なので、ζ は Γ の全元で固定される
  have hfix : ∀ g : (selfField p).absGal, g ζ = ζ := by
    intro g
    have hs := cyclotomicCharacter.spec (L := (selfField p).closure) p g.toRingEquiv ζ hζpow
    rw [hcon g] at hs
    simpa [ZMod.val_one] using hs
  -- 3. `K̄/ℚ_p` は Galois なので `fixedField ⊤ = ⊥`、すなわち ζ ∈ ℚ_p
  have hmem : ζ ∈ Set.range (algebraMap (selfField p).carrier (selfField p).closure) := by
    rw [← IntermediateField.mem_bot,
      ← InfiniteGalois.fixedField_fixingSubgroup
        (⊥ : IntermediateField (selfField p).carrier (selfField p).closure),
      IntermediateField.fixingSubgroup_bot, IntermediateField.mem_fixedField_iff]
    intro f _
    exact hfix f
  obtain ⟨y, hy⟩ := hmem
  have hinj := (algebraMap (selfField p).carrier (selfField p).closure).injective
  have hy2 : y ^ p ^ 2 = 1 := by
    apply hinj; rw [map_pow, hy, map_one]; exact hζpow
  -- 4. これは偽
  rcases Nat.lt_or_ge p 3 with hlt | hp3
  · -- `p = 2`: ζ は原始 4 乗根なので `ζ² ≠ 1`、しかし `ℚ_2` では `ζ² = 1`
    have hp2 : p = 2 := by have := hp.out.two_le; omega
    subst hp2
    have h4 : y ^ 4 = 1 := by rw [show (4 : ℕ) = 2 ^ 2 from rfl]; exact hy2
    have hy1 : y ^ 2 = 1 := sq_eq_one_of_pow_four_eq_one_q2 h4
    refine hζ.pow_ne_one_of_pos_of_lt (l := 2) (by norm_num) (by norm_num) ?_
    rw [← hy, ← map_pow, hy1, map_one]
  · -- `p ≥ 3`: `y^p` は `ℚ_p` の中の `p` 乗根なので 1、しかし `ζ^p ≠ 1`
    have hηp : (y ^ p) ^ p = 1 := by rw [← pow_mul, ← pow_two]; exact hy2
    have hη : y ^ p = 1 := eq_one_of_pow_p_eq_one_qp p hp3 hηp
    refine hζ.pow_ne_one_of_pos_of_lt (l := p) (by omega) (by nlinarith) ?_
    rw [← hy, ← map_pow, hη, map_one]

/-- **`cyclotomicCharacterObject` は定数関数 1 を返さない。**

Proposition 1.1(`cyclotomicCharacter_recoverable`)が主張しているのは
「`transport α (obj K) = obj K'`」であって、もし `obj` が常に定数 1 なら
両辺が定数 1 で自明に成り立つ。この 1 行がその退化を塞ぐ。 -/
theorem cyclotomicCharacterObject_selfField_ne_const_one (p : ℕ) [Fact p.Prime] :
    (cyclotomicCharacterObject (p := p)).obj (selfField p) ≠ fun _ => 1 := by
  intro hcon
  obtain ⟨g, hg⟩ := cyclotomicCharacter_selfField_ne_one p
  exact hg (congrFun hcon g)

end ABC3.Check.PGC
