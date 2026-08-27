import ABC3.Found.GenEll.Conductor
import Mathlib.NumberTheory.RamificationInertia.Basic

/-!
# [GenEll] Remark 1.5.1 / Proposition 1.7 —— 「`Σ` の上の寄与は `≈ 0`」(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9–p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★★★★これが BD-class を使う理由そのものである

`Remark 1.5.1` は「別の `ℤ`-モデルへの同型が**ある有限素数集合 `Σ` の外で**延びる」ことから
`log-cond` の BD-class が `(X_ℚ, D_ℚ)` だけに依る、と言う。
`Proposition 1.7` の証明も『`Σ` の上の寄与は `≈ 0`』として同じ形を使う。

★**その「`≈ 0`」の中身がこれである**:

    log-cond の Σ の上の部分 ≤ Σ_{q ∈ Σ} log q

★★右辺は **`Σ` だけで決まる定数**であり、**点 `x` にも定義体 `F` にも依らない**。
だから差が有界(= BD-同値)になる。

## ★★★★一様性はどこから来るか

`log-cond` は `deg` を `[F:ℚ]` で割った**正規化次数**である(`degNormalized`)。
`q` の上の素点 `v` について `N v = q^{f_v}` なので、`Σ` の上の寄与は

    (Σ_{v | q} f_v) · log q / [F:ℚ]

★★分子の `Σ_{v|q} f_v` は `Σ_{v|q} e_v f_v = [F:ℚ]` で抑えられる(`e_v ≥ 1`)。
**ここで `[F:ℚ]` が約分される**——これが一様性の機構である。

★★★基本等式 `Σ e f = [L:K]` は mathlib にある
(`Ideal.sum_ramification_inertia`)。★本ファイルはそこから
`Σ_{v|q} f_v ≤ [F:ℚ]` を取り出して `log` の側へ運ぶ。
-/

namespace ABC3.Found.GenEll

open scoped BigOperators

variable {F : Type*} [Field F] [NumberField F]

/-! ## ★被約算術因子の次数は台の上の `log N v` の和 -/

/-- ★**被約因子の次数は台の上で `log N v` を足したもの**(係数が `0` か `1` だから)。 -/
theorem deg_adivRed_le_sum (a : ADiv F) (V : Finset (FinitePlace F))
    (hsupp : (ADivRed a).fin.support ⊆ V) :
    deg (ADivRed a) ≤ ∑ v ∈ V, Real.log (residueCard v) := by
  classical
  have harc : (ADivRed a).arc.sum (fun _ r => r) = 0 := by
    rw [ADivRed_arc]; simp
  rw [deg, harc, add_zero, Finsupp.sum]
  have hterm : ∀ v ∈ (ADivRed a).fin.support,
      (((ADivRed a).fin v : ℤ) : ℝ) * Real.log (residueCard v) = Real.log (residueCard v) := by
    intro v hv
    rw [Finsupp.mem_support_iff, ADivRed_fin] at hv
    have h1 : (ADivRed a).fin v = 1 := by
      rw [ADivRed_fin]
      split_ifs at hv ⊢ with h
      · rfl
      · exact absurd rfl hv
    rw [h1]
    push_cast
    ring
  rw [Finset.sum_congr rfl hterm]
  refine Finset.sum_le_sum_of_subset_of_nonneg hsupp (fun v _ _ => ?_)
  exact (log_residueCard_pos v).le

/-! ## ★★★★★★`Σ_{v|q} f_v ≤ [F:ℚ]` —— 基本等式から -/

/-- ★★★★★★**惰性次数の和は拡大次数を超えない**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★`Σ_{v|p} e_v f_v = [F:ℚ]`(`Ideal.sum_ramification_inertia`)と `e_v ≥ 1` から。
★★**これが一様性の源**である——右辺が `[F:ℚ]` なので、
正規化(`/[F:ℚ]`)で消える。 -/
theorem sum_inertiaDeg_le (p : Ideal ℤ) [p.IsMaximal] (hp0 : p ≠ ⊥)
    (V : Finset (FinitePlace F))
    (hV : ∀ v ∈ V, (v.asIdeal).LiesOver p) :
    ∑ v ∈ V, p.inertiaDeg v.asIdeal ≤ Module.finrank ℚ F := by
  classical
  have hinj : Set.InjOn (fun v : FinitePlace F => v.asIdeal) ↑V := by
    intro a _ b _ h
    exact IsDedekindDomain.HeightOneSpectrum.ext h
  have himg : V.image (fun v : FinitePlace F => v.asIdeal)
      ⊆ IsDedekindDomain.primesOverFinset p (NumberField.RingOfIntegers F) := by
    intro P hP
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.1 hP
    exact (IsDedekindDomain.mem_primesOverFinset_iff hp0 _).2 ⟨v.isPrime, hV v hv⟩
  calc ∑ v ∈ V, p.inertiaDeg v.asIdeal
      = ∑ P ∈ V.image (fun v : FinitePlace F => v.asIdeal), p.inertiaDeg P :=
        (Finset.sum_image hinj).symm
    _ ≤ ∑ P ∈ IsDedekindDomain.primesOverFinset p (NumberField.RingOfIntegers F),
          p.inertiaDeg P := Finset.sum_le_sum_of_subset himg
    _ ≤ ∑ P ∈ IsDedekindDomain.primesOverFinset p (NumberField.RingOfIntegers F),
          p.ramificationIdx P * p.inertiaDeg P := by
        refine Finset.sum_le_sum (fun P hP => ?_)
        have h1 : p.ramificationIdx P ≠ 0 := Ideal.Factors.ramificationIdx_ne_zero p ⟨P, hP⟩
        exact Nat.le_mul_of_pos_left _ (Nat.pos_of_ne_zero h1)
    _ = Module.finrank ℚ F := Ideal.sum_ramification_inertia _ ℚ F hp0

/-! ## ★★★`N v = q^{f_v}` を `log` の側へ -/

/-- ★★**`log N v = f_v · log q`**(`v` が `q` の上にあるとき)。 -/
theorem log_residueCard_eq_inertiaDeg_mul (q : ℕ) (hq : q.Prime) (v : FinitePlace F)
    [v.asIdeal.LiesOver (Ideal.span {(q : ℤ)})] :
    Real.log (residueCard v)
      = ((Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal : ℝ) * Real.log q := by
  have hnorm : residueCard v = q ^ (Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal := by
    have := Ideal.absNorm_eq_pow_inertiaDeg (R := NumberField.RingOfIntegers F)
      (p := (q : ℤ)) v.asIdeal (Nat.prime_iff_prime_int.mp hq)
    simpa [residueCard, Int.natAbs_natCast] using this
  rw [hnorm]
  push_cast
  rw [Real.log_pow]

/-- ★★★★★**1 つの素数 `q` の上の寄与は `[F:ℚ]·log q` 以下**。 -/
theorem sum_log_residueCard_le (q : ℕ) (hq : q.Prime) (V : Finset (FinitePlace F))
    (hV : ∀ v ∈ V, (v.asIdeal).LiesOver (Ideal.span {(q : ℤ)})) :
    ∑ v ∈ V, Real.log (residueCard v) ≤ (Module.finrank ℚ F : ℝ) * Real.log q := by
  classical
  have h0 : (q : ℤ) ≠ 0 := by exact_mod_cast hq.ne_zero
  have hp0 : (Ideal.span {(q : ℤ)} : Ideal ℤ) ≠ ⊥ := by simpa using h0
  haveI hmax : (Ideal.span {(q : ℤ)}).IsMaximal :=
    ((Ideal.span_singleton_prime h0).mpr (Nat.prime_iff_prime_int.mp hq)).isMaximal hp0
  have hsum : ∑ v ∈ V, Real.log (residueCard v)
      = (∑ v ∈ V, ((Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal : ℕ) : ℝ) * Real.log q := by
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl (fun v hv => ?_)
    haveI := hV v hv
    rw [log_residueCard_eq_inertiaDeg_mul q hq v]
  rw [hsum]
  have hle : (∑ v ∈ V, (((Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal : ℕ) : ℝ))
      ≤ (Module.finrank ℚ F : ℝ) := by
    have hnat := sum_inertiaDeg_le (Ideal.span {(q : ℤ)}) hp0 V hV
    calc (∑ v ∈ V, (((Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal : ℕ) : ℝ))
        = ((∑ v ∈ V, ((Ideal.span {(q : ℤ)}).inertiaDeg v.asIdeal) : ℕ) : ℝ) := by push_cast; ring
      _ ≤ (Module.finrank ℚ F : ℝ) := by exact_mod_cast hnat
  have hlog : 0 ≤ Real.log q := Real.log_nonneg (by exact_mod_cast hq.one_lt.le)
  exact mul_le_mul_of_nonneg_right hle hlog

/-! ## ★★★★★★★★★「`Σ` の上の寄与は `≈ 0`」 -/

/-- ★★★★★★★★★**`Σ` の上の `log-cond` の寄与は `Σ_{q ∈ Σ} log q` で抑えられる**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★★**右辺は `Σ` だけで決まる定数**であり、点にも定義体 `F` にも依らない
——`[F:ℚ]` は正規化で約分される。★★★これが BD-class が吸収する「緩み」である。 -/
theorem degNormalized_adivRed_le_sum_log (a : ADiv F) (Sig : Finset ℕ)
    (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ (ADivRed a).fin.support, ch v ∈ Sig)
    (hover : ∀ v ∈ (ADivRed a).fin.support,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    degNormalized (ADivRed a) ≤ ∑ q ∈ Sig, Real.log q := by
  classical
  have hn : 0 < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hdeg : deg (ADivRed a)
      ≤ ∑ v ∈ (ADivRed a).fin.support, Real.log (residueCard v) :=
    deg_adivRed_le_sum a _ (Finset.Subset.refl _)
  have hfib : ∑ q ∈ Sig, ∑ v ∈ (ADivRed a).fin.support with ch v = q,
        Real.log (residueCard v)
      = ∑ v ∈ (ADivRed a).fin.support, Real.log (residueCard v) :=
    Finset.sum_fiberwise_of_maps_to hmem _
  have hq : ∀ q ∈ Sig, ∑ v ∈ (ADivRed a).fin.support with ch v = q,
      Real.log (residueCard v) ≤ (Module.finrank ℚ F : ℝ) * Real.log q := by
    intro q hqS
    refine sum_log_residueCard_le q (hprime q hqS) _ (fun v hv => ?_)
    obtain ⟨hv1, hv2⟩ := Finset.mem_filter.1 hv
    have := hover v hv1
    rwa [hv2] at this
  have hsum : ∑ v ∈ (ADivRed a).fin.support, Real.log (residueCard v)
      ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := by
    rw [← hfib]
    exact Finset.sum_le_sum hq
  rw [degNormalized, div_le_iff₀ hn]
  calc deg (ADivRed a) ≤ ∑ v ∈ (ADivRed a).fin.support, Real.log (residueCard v) := hdeg
    _ ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := hsum
    _ = (∑ q ∈ Sig, Real.log q) * (Module.finrank ℚ F : ℝ) := by
        rw [Finset.sum_mul]
        exact Finset.sum_congr rfl (fun q _ => by ring)

/-! ## ★★★★★★★★★★一般形 —— 係数が `≤ 1` の因子なら何でもよい

★`ADivRed` の形に限る必要はない。`log-cond` の導手 `f_x^D ≝ (D_x)_red` は
`CartierPullback.lean` では **イデアルの根基**として作られており、
`ADivRed` を経由していない。★★一般形にしておけば、そちらへも直接使える。
-/

/-- ★★**係数が `≤ 1` の因子の次数は台の上の `log N v` の和で抑えられる**。 -/
theorem deg_le_sum_of_coeff_le_one (c : ADiv F) (harc : c.arc = 0)
    (hle1 : ∀ v, c.fin v ≤ 1)
    (V : Finset (FinitePlace F)) (hsupp : c.fin.support ⊆ V) :
    deg c ≤ ∑ v ∈ V, Real.log (residueCard v) := by
  classical
  have h0 : c.arc.sum (fun _ r => r) = 0 := by rw [harc]; simp
  rw [deg, h0, add_zero, Finsupp.sum]
  have hterm : ∀ v ∈ c.fin.support,
      ((c.fin v : ℤ) : ℝ) * Real.log (residueCard v) ≤ Real.log (residueCard v) := by
    intro v _
    have h1 : ((c.fin v : ℤ) : ℝ) ≤ 1 := by exact_mod_cast hle1 v
    nlinarith [(log_residueCard_pos v).le, h1]
  refine le_trans (Finset.sum_le_sum hterm) ?_
  exact Finset.sum_le_sum_of_subset_of_nonneg hsupp
    (fun v _ _ => (log_residueCard_pos v).le)

/-- ★★★★★★★★★★**一般形の `Σ`-限界**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

係数が `≤ 1`、台が `Σ` の上にあるなら、正規化次数は `Σ_{q ∈ Σ} log q` 以下。
★**`ADivRed` を経由しない導手(イデアルの根基)にもそのまま使える**。 -/
theorem degNormalized_le_sum_log (c : ADiv F) (harc : c.arc = 0)
    (hle1 : ∀ v, c.fin v ≤ 1)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ c.fin.support, ch v ∈ Sig)
    (hover : ∀ v ∈ c.fin.support, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    degNormalized c ≤ ∑ q ∈ Sig, Real.log q := by
  classical
  have hn : 0 < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hdeg : deg c ≤ ∑ v ∈ c.fin.support, Real.log (residueCard v) :=
    deg_le_sum_of_coeff_le_one c harc hle1 _ (Finset.Subset.refl _)
  have hfib : ∑ q ∈ Sig, ∑ v ∈ c.fin.support with ch v = q, Real.log (residueCard v)
      = ∑ v ∈ c.fin.support, Real.log (residueCard v) :=
    Finset.sum_fiberwise_of_maps_to hmem _
  have hq : ∀ q ∈ Sig, ∑ v ∈ c.fin.support with ch v = q,
      Real.log (residueCard v) ≤ (Module.finrank ℚ F : ℝ) * Real.log q := by
    intro q hqS
    refine sum_log_residueCard_le q (hprime q hqS) _ (fun v hv => ?_)
    obtain ⟨hv1, hv2⟩ := Finset.mem_filter.1 hv
    have hh := hover v hv1
    rwa [hv2] at hh
  have hsum : ∑ v ∈ c.fin.support, Real.log (residueCard v)
      ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := by
    rw [← hfib]
    exact Finset.sum_le_sum hq
  rw [degNormalized, div_le_iff₀ hn]
  calc deg c ≤ ∑ v ∈ c.fin.support, Real.log (residueCard v) := hdeg
    _ ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := hsum
    _ = (∑ q ∈ Sig, Real.log q) * (Module.finrank ℚ F : ℝ) := by
        rw [Finset.sum_mul]
        exact Finset.sum_congr rfl (fun q _ => by ring)

/-! ## ★★★★★★★★★★★2 つの因子の差 —— BD-class の中身 -/

/-- ★★**台が `Σ` の上にある集まりの `log N v` の総和の限界**(上の証明の中核を切り出したもの)。 -/
theorem sum_log_residueCard_le_of_over (V : Finset (FinitePlace F)) (Sig : Finset ℕ)
    (hprime : ∀ q ∈ Sig, q.Prime) (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ V, ch v ∈ Sig)
    (hover : ∀ v ∈ V, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    ∑ v ∈ V, Real.log (residueCard v)
      ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
  classical
  have hfib : ∑ q ∈ Sig, ∑ v ∈ V with ch v = q, Real.log (residueCard v)
      = ∑ v ∈ V, Real.log (residueCard v) :=
    Finset.sum_fiberwise_of_maps_to hmem _
  have hq : ∀ q ∈ Sig, ∑ v ∈ V with ch v = q,
      Real.log (residueCard v) ≤ (Module.finrank ℚ F : ℝ) * Real.log q := by
    intro q hqS
    refine sum_log_residueCard_le q (hprime q hqS) _ (fun v hv => ?_)
    obtain ⟨hv1, hv2⟩ := Finset.mem_filter.1 hv
    have hh := hover v hv1
    rwa [hv2] at hh
  calc ∑ v ∈ V, Real.log (residueCard v)
      = ∑ q ∈ Sig, ∑ v ∈ V with ch v = q, Real.log (residueCard v) := hfib.symm
    _ ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := Finset.sum_le_sum hq
    _ = (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by rw [Finset.mul_sum]

/-- ★★★★★★★★★★★**`Σ` の外で一致する 2 つの因子の次数の差は `Σ` だけで決まる定数以下**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX(Q) depends only on the pair

★★これが `Remark 1.5.1` の『BD-class は `(X_ℚ, D_ℚ)` だけに依る』の**計算部分**である
——`Σ` の外で導手が一致すれば、正規化次数の差は `Σ_{q∈Σ} log q` を超えない。
★★★**右辺は点にも定義体にも依らない**。だから差が有界(= BD-同値)になる。 -/
theorem degNormalized_sub_le_sum_log (c c' : ADiv F)
    (harc : c.arc = 0) (harc' : c'.arc = 0)
    (hc1 : ∀ v, c.fin v ≤ 1) (hc'0 : ∀ v, 0 ≤ c'.fin v)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hagree : ∀ v, ch v ∉ Sig → c.fin v = c'.fin v)
    (hover : ∀ v : FinitePlace F, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    degNormalized c - degNormalized c' ≤ ∑ q ∈ Sig, Real.log q := by
  classical
  have hn : 0 < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  set T : Finset (FinitePlace F) := c.fin.support ∪ c'.fin.support with hT
  have hgz : ∀ (v : FinitePlace F), ((0 : ℤ) : ℝ) * Real.log (residueCard v) = 0 := by
    intro v; simp
  have hdegc : deg c = ∑ v ∈ T, ((c.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    rw [deg, harc]
    simp only [Finsupp.sum_zero_index, add_zero]
    exact Finsupp.sum_of_support_subset _ Finset.subset_union_left _ (fun v _ => hgz v)
  have hdegc' : deg c' = ∑ v ∈ T, ((c'.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    rw [deg, harc']
    simp only [Finsupp.sum_zero_index, add_zero]
    exact Finsupp.sum_of_support_subset _ Finset.subset_union_right _ (fun v _ => hgz v)
  have hdiff : deg c - deg c'
      = ∑ v ∈ T, (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v) := by
    rw [hdegc, hdegc', ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun v _ => by ring)
  -- `Σ` の外の項は消える
  have hrestrict : ∑ v ∈ T,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v)
      = ∑ v ∈ T with ch v ∈ Sig,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v) := by
    refine (Finset.sum_subset (Finset.filter_subset _ _) (fun v hv hnot => ?_)).symm
    have hns : ch v ∉ Sig := by
      intro hc
      exact hnot (Finset.mem_filter.2 ⟨hv, hc⟩)
    rw [hagree v hns]
    ring
  -- 各項を `log N v` で抑える
  have hbound : ∑ v ∈ T with ch v ∈ Sig,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v)
      ≤ ∑ v ∈ T with ch v ∈ Sig, Real.log (residueCard v) := by
    refine Finset.sum_le_sum (fun v _ => ?_)
    have h1 : ((c.fin v : ℤ) : ℝ) ≤ 1 := by exact_mod_cast hc1 v
    have h2 : (0 : ℝ) ≤ ((c'.fin v : ℤ) : ℝ) := by exact_mod_cast hc'0 v
    nlinarith [(log_residueCard_pos v).le]
  have hSig : ∑ v ∈ T with ch v ∈ Sig, Real.log (residueCard v)
      ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
    refine sum_log_residueCard_le_of_over _ Sig hprime ch (fun v hv => (Finset.mem_filter.1 hv).2)
      (fun v _ => hover v)
  have hkey : deg c - deg c' ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
    rw [hdiff, hrestrict]
    exact le_trans hbound hSig
  rw [degNormalized, degNormalized, div_sub_div_same, div_le_iff₀ hn]
  linarith [hkey]

/-! ## ★出典の紐付け(`.src`) -/

def sum_inertiaDeg_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の寄与が ≈ 0 であること——惰性次数の和の限界)",
    sectionId := "genell-prop-1-7" }

def degNormalized_adivRed_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の log-cond の寄与は Σ_{q∈Σ} log q で抑えられる)",
    sectionId := "genell-prop-1-7" }

def degNormalized_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の寄与の限界の一般形——係数が ≤ 1 なら何でもよい)",
    sectionId := "genell-prop-1-7" }

def degNormalized_sub_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(BD-class が (X_ℚ, D_ℚ) だけに依ることの計算部分)",
    sectionId := "genell-rem-1-5-1" }

end ABC3.Found.GenEll
