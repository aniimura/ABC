import ABC3.Found.PGC.LubinTateFormalGroupLaw

/-!
# Lubin-Tate 形式群法則: 一意性補題への準備(進行中)

古典的な Lubin-Tate の一意性補題(2つの1変数冪級数が同じ次数1の係数を
持ち、同じ関数等式 `φ(f(X))=f(φ(X))` を満たせば等しい)を経由して
`F_f(X,0)=X`(単位元則)を示す計画の、最初の部品——1変数冪級数を
「`X_1` に依存しない2変数冪級数」として埋め込む写像 `emb` の基本性質。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A]

/-- 1変数冪級数を `X_1` に依存しない2変数冪級数として埋め込む
(`T ↦ X_0`)。 -/
noncomputable def emb (p : PowerSeries A) : MvPowerSeries (Fin 2) A :=
  PowerSeries.subst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) p

theorem hasSubst_X0 : PowerSeries.HasSubst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) := by
  show IsNilpotent (MvPowerSeries.constantCoeff (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A))
  rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero

theorem emb_add (p q : PowerSeries A) : emb (p + q) = emb p + emb q :=
  PowerSeries.subst_add hasSubst_X0 p q

theorem emb_sub (p q : PowerSeries A) : emb (p - q) = emb p - emb q :=
  PowerSeries.subst_sub hasSubst_X0 p q

/-- `emb p` の係数——`X_0` の次数 `n`・`X_1` の次数 `0` の単項式の係数は
`p` の `n` 次係数、それ以外は `0`。 -/
theorem coeff_emb (p : PowerSeries A) (e : Fin 2 →₀ ℕ) :
    MvPowerSeries.coeff e (emb p) =
      if e 1 = 0 then PowerSeries.coeff (e 0) p else 0 := by
  rw [emb, PowerSeries.coeff_subst hasSubst_X0]
  rw [finsum_eq_sum_of_support_subset _ (s := {e 0}) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    apply hd
    have : MvPowerSeries.coeff e ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) ^ d) = 0 := by
      rw [MvPowerSeries.coeff_X_pow]
      split_ifs with hcase
      · exfalso; apply hcon
        rw [hcase, Finsupp.single_eq_same]
      · rfl
    simp [this])]
  by_cases h1 : e 1 = 0
  · rw [if_pos h1]
    have he : e = Finsupp.single (0 : Fin 2) (e 0) := by
      ext i
      fin_cases i
      · simp
      · simp [h1]
    rw [Finset.sum_singleton, he, MvPowerSeries.coeff_X_pow]
    simp
  · rw [if_neg h1]
    rw [Finset.sum_singleton]
    have : MvPowerSeries.coeff e ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) ^ (e 0)) = 0 := by
      rw [MvPowerSeries.coeff_X_pow]
      split_ifs with hcase
      · exfalso; apply h1
        rw [hcase, Finsupp.single_eq_of_ne' (by decide : (0 : Fin 2) ≠ 1)]
      · rfl
    simp [this]

/-- `p` の次数が `≥n` なら `emb p` の次数も `≥n`。 -/
theorem order_emb_ge {p : PowerSeries A} {n : ℕ} (hp : (n : ℕ∞) ≤ p.order) :
    (n : ℕ∞) ≤ (emb p).order := by
  apply MvPowerSeries.nat_le_order
  intro d hd
  rw [coeff_emb]
  by_cases h1 : d 1 = 0
  · rw [if_pos h1]
    apply PowerSeries.coeff_of_lt_order
    have hdeg : Finsupp.degree d = d 0 := by
      rw [finsupp_degree_fin2, h1, add_zero]
    rw [hdeg] at hd
    calc ((d 0 : ℕ) : ℕ∞) < (n : ℕ∞) := by exact_mod_cast hd
      _ ≤ p.order := hp
  · rw [if_neg h1]

theorem emb_injective : Function.Injective (emb (A := A)) := by
  intro p q hpq
  ext n
  have := coeff_emb (A := A) p (Finsupp.single (0 : Fin 2) n)
  have hq := coeff_emb (A := A) q (Finsupp.single (0 : Fin 2) n)
  simp only [Finsupp.single_eq_same] at this hq
  rw [hpq] at this
  rw [this] at hq
  simpa using hq

/-- `emb` は代入と両立する: `emb(p)` を `f` に代入した結果は、`p` を `f` に
代入してから埋め込んだ結果に一致する——`emb(f(p(T))) = f(emb(p))`。
`PowerSeries.subst_comp_subst_apply` を `a:=p`・`b:=X_0` で適用するだけ。 -/
theorem subst_emb (p f : PowerSeries A) (hp : PowerSeries.HasSubst p) :
    PowerSeries.subst (emb p) f = emb (PowerSeries.subst p f) :=
  (PowerSeries.subst_comp_subst_apply hp hasSubst_X0 f).symm

/-- `h ≡ πX (mod deg 2)` のとき、`h^{n+1}` の次数 `n+1` の係数はちょうど
`π^{n+1}`(先頭項どうしの積だけが効く)。 -/
theorem coeff_pow_eq_of_order_one {h : PowerSeries A} {π : A}
    (hh0 : PowerSeries.constantCoeff h = 0) (hh1 : PowerSeries.coeff 1 h = π) (n : ℕ) :
    PowerSeries.coeff (n + 1) (h ^ (n + 1)) = π ^ (n + 1) := by
  induction n with
  | zero => simpa using hh1
  | succ n ih =>
    have horder : ((n + 1 : ℕ) : ℕ∞) ≤ (h ^ (n + 1)).order := by
      calc ((n + 1 : ℕ) : ℕ∞) = (n + 1) • (1 : ℕ∞) := by simp
        _ ≤ (n + 1) • h.order := by
            gcongr
            exact PowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hh0
        _ ≤ (h ^ (n + 1)).order := PowerSeries.le_order_pow h (n + 1)
    have hstep : h ^ (n + 1 + 1) = h ^ (n + 1) * h := by ring
    rw [hstep, PowerSeries.coeff_mul]
    rw [Finset.sum_eq_single (n + 1, 1)]
    · rw [ih, hh1]; ring
    · intro b hb hbne
      simp only [Finset.mem_antidiagonal] at hb
      by_cases hb1 : b.1 < n + 1
      · have : PowerSeries.coeff b.1 (h ^ (n + 1)) = 0 := by
          apply PowerSeries.coeff_of_lt_order
          exact lt_of_lt_of_le (by exact_mod_cast hb1) horder
        simp [this]
      · have hb2le : b.2 ≤ 1 := by omega
        interval_cases hb2v : b.2
        · simp [hh0]
        · exfalso; apply hbne
          have hb1eq : b.1 = n + 1 := by omega
          exact Prod.ext hb1eq hb2v
    · intro hnotmem
      exfalso; apply hnotmem
      simp [Finset.mem_antidiagonal]

/-- `δ` の次数が `≥n+1` のとき、`δ` を `h`(`≡πX(mod deg2)`)に代入した結果の
次数 `n+1` の係数は `coeff(n+1)δ・π^{n+1}`——`k=n+1` の項だけが効く。 -/
theorem coeff_subst_eq_of_order_ge {δ h : PowerSeries A} {π : A} {n : ℕ}
    (hh0 : PowerSeries.constantCoeff h = 0) (hh1 : PowerSeries.coeff 1 h = π)
    (hδorder : ((n + 1 : ℕ) : ℕ∞) ≤ δ.order) (hHS : PowerSeries.HasSubst h) :
    PowerSeries.coeff (n + 1) (PowerSeries.subst h δ) =
      PowerSeries.coeff (n + 1) δ * π ^ (n + 1) := by
  rw [PowerSeries.coeff_subst' hHS]
  rw [finsum_eq_sum_of_support_subset _ (s := ({n + 1} : Finset ℕ)) (fun k hk => by
    simp only [Function.mem_support] at hk
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    apply hk
    rcases lt_or_gt_of_ne hcon with hlt | hgt
    · have : PowerSeries.coeff k δ = 0 := by
        apply PowerSeries.coeff_of_lt_order
        exact lt_of_lt_of_le (by exact_mod_cast hlt) hδorder
      simp [this]
    · have horderh : (1 : ℕ∞) ≤ h.order :=
        PowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hh0
      have hpow : (k : ℕ∞) ≤ (h ^ k).order := by
        calc (k : ℕ∞) = k • (1 : ℕ∞) := by simp
          _ ≤ k • h.order := by gcongr
          _ ≤ (h ^ k).order := PowerSeries.le_order_pow h k
      have : PowerSeries.coeff (n + 1) (h ^ k) = 0 := by
        apply PowerSeries.coeff_of_lt_order
        calc ((n + 1 : ℕ) : ℕ∞) < (k : ℕ∞) := by exact_mod_cast hgt
          _ ≤ _ := hpow
      simp [this])]
  simp [coeff_pow_eq_of_order_one hh0 hh1 n]

/-- ★★★★★★★★★**Lubin-Tate の一意性補題**。`h≡πX(mod deg2)` に対し、次数1の
係数が一致し(かつ定数項0)、同じ関数等式 `h(α(T))=α(h(T))`・`h(β(T))=β(h(T))`
を満たす2つの冪級数 `α,β` は一致する——古典的な Lubin-Tate 理論の基礎
(Lubin–Tate 1965)。`α−β` の次数ごとの帰納法で示す: `emb`・`subst_emb`・
`coeff_subst_linearize`(2変数の線形化)・`coeff_subst_eq_of_order_ge`
(合成の先頭項)を組み合わせ、`coeff(n+1)(α−β)・(π^{n+1}−π)=0` を得て、
`π^n≠1`(`π` が極大イデアルの元、`n≥1`)から整域で割り切る。 -/
theorem powerSeries_uniqueness [IsLocalRing A] [IsDomain A] {π : A}
    (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {h : PowerSeries A} (hh0 : PowerSeries.constantCoeff h = 0) (hh1 : PowerSeries.coeff 1 h = π)
    {α β : PowerSeries A} (hα0 : PowerSeries.constantCoeff α = 0)
    (hβ0 : PowerSeries.constantCoeff β = 0)
    (hlead : PowerSeries.coeff 1 α = PowerSeries.coeff 1 β)
    (heqα : PowerSeries.subst h α = PowerSeries.subst α h)
    (heqβ : PowerSeries.subst h β = PowerSeries.subst β h) :
    α = β := by
  have hHSh : PowerSeries.HasSubst h := by
    show IsNilpotent (PowerSeries.constantCoeff h); rw [hh0]; exact IsNilpotent.zero
  have hHSα : PowerSeries.HasSubst α := by
    show IsNilpotent (PowerSeries.constantCoeff α); rw [hα0]; exact IsNilpotent.zero
  have hHSβ : PowerSeries.HasSubst β := by
    show IsNilpotent (PowerSeries.constantCoeff β); rw [hβ0]; exact IsNilpotent.zero
  set δ := α - β with hδ_def
  have hδ0 : PowerSeries.constantCoeff δ = 0 := by rw [hδ_def, map_sub, hα0, hβ0, sub_zero]
  have hδ1 : PowerSeries.coeff 1 δ = 0 := by
    rw [hδ_def, map_sub, hlead, sub_self]
  have hbase : ((2 : ℕ) : ℕ∞) ≤ δ.order := by
    apply PowerSeries.nat_le_order
    intro i hi
    interval_cases i
    · rw [PowerSeries.coeff_zero_eq_constantCoeff]; exact hδ0
    · exact hδ1
  have hstep : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order →
      ((n + 2 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn1 hδorder
    apply PowerSeries.nat_le_order
    intro i hi
    rcases lt_or_eq_of_le (by omega : i ≤ n + 1) with hilt | hieq
    · exact PowerSeries.coeff_of_lt_order i (lt_of_lt_of_le (by exact_mod_cast hilt) hδorder)
    · subst hieq
      set e : Fin 2 →₀ ℕ := Finsupp.single (0 : Fin 2) (n + 1) with he_def
      have he0 : (e 0 : ℕ) = n + 1 := by rw [he_def, Finsupp.single_eq_same]
      have he1 : (e 1 : ℕ) = 0 := by
        rw [he_def, Finsupp.single_eq_of_ne' (by decide : (0 : Fin 2) ≠ 1)]
      have hedeg : Finsupp.degree e = n + 1 := by rw [finsupp_degree_fin2, he0, he1, add_zero]
      have hh0' : PowerSeries.coeff 0 h = 0 := by
        rw [PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hh0
      have hΦ0 : MvPowerSeries.constantCoeff (emb β) = 0 := by
        rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_emb]
        simp [PowerSeries.coeff_zero_eq_constantCoeff_apply, hβ0]
      have hφ0 : MvPowerSeries.constantCoeff (emb δ) = 0 := by
        rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_emb]
        simp [PowerSeries.coeff_zero_eq_constantCoeff_apply, hδ0]
      have hΦorder : (1 : ℕ∞) ≤ (emb β).order :=
        MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hΦ0
      have hφorder : ((n + 1 : ℕ) : ℕ∞) ≤ (emb δ).order := order_emb_ge hδorder
      have hlin := coeff_subst_linearize hΦ0 hφ0 hΦorder hφorder (by omega : 1 ≤ n + 1) h π hh0' hh1 e
        (by rw [hedeg])
      have hsum : emb β + emb δ = emb α := by
        rw [← emb_add]; congr 1; rw [hδ_def]; ring
      rw [hsum] at hlin
      rw [subst_emb α h hHSα, subst_emb β h hHSβ] at hlin
      rw [← map_sub (MvPowerSeries.coeff e), ← emb_sub] at hlin
      rw [show PowerSeries.subst α h - PowerSeries.subst β h =
          PowerSeries.subst h α - PowerSeries.subst h β from by rw [heqα, heqβ]] at hlin
      rw [← PowerSeries.subst_sub hHSh α β] at hlin
      have hcoeffe : ∀ p : PowerSeries A, MvPowerSeries.coeff e (emb p) = PowerSeries.coeff (n + 1) p := by
        intro p
        rw [coeff_emb, he0, he1]
        simp
      rw [hcoeffe, hcoeffe] at hlin
      have hkey := coeff_subst_eq_of_order_ge (δ := δ) (h := h) (π := π) (n := n) hh0 hh1 hδorder hHSh
      rw [hkey] at hlin
      have hfactor : PowerSeries.coeff (n + 1) δ * (π ^ (n + 1) - π) = 0 := by
        have := hlin
        ring_nf
        ring_nf at this
        linear_combination this
      have hπn_mem : π ^ n ∈ IsLocalRing.maximalIdeal A :=
        Ideal.pow_mem_of_mem _ (hπmax ▸ Ideal.mem_span_singleton_self π) n (by omega)
      have hπn_ne_one : π ^ n ≠ 1 := fun heq => by
        rw [heq] at hπn_mem
        exact IsLocalRing.maximalIdeal.isMaximal A |>.ne_top
          (Ideal.eq_top_of_isUnit_mem _ hπn_mem isUnit_one)
      have hne : π ^ (n + 1) - π ≠ 0 := by
        intro hcon
        apply hπn_ne_one
        have : π * (π ^ n - 1) = 0 := by ring_nf; ring_nf at hcon; linear_combination hcon
        rcases mul_eq_zero.mp this with h1 | h2
        · exact absurd h1 hπne0
        · exact sub_eq_zero.mp h2
      exact (mul_eq_zero.mp hfactor).resolve_right hne
  have hall : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn
    induction n, hn using Nat.le_induction with
    | base => exact hbase
    | succ n hn1 ih => exact hstep n hn1 ih
  have horder : δ.order = ⊤ := by
    by_contra hne
    obtain ⟨m, hm⟩ := WithTop.ne_top_iff_exists.mp hne
    have := hall (m + 1) (by omega)
    rw [← hm] at this
    have hcontra : ((m : ℕ) : ℕ∞) < ((m + 1 + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega)
    exact absurd this (not_le.mpr hcontra)
  have hδzero : δ = 0 := PowerSeries.order_eq_top.mp horder
  have := hδzero
  rw [hδ_def, sub_eq_zero] at this
  exact this

end ABC3.Found.PGC
