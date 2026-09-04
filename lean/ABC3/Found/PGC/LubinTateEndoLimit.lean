import ABC3.Found.PGC.LubinTateEndoSequence

/-!
# `𝒪_K` 作用への拡張(1変数版)の近似列の極限(`sorry` 無し)

`Found/PGC/LubinTateEndoSequence.lean::φSeq_endo` の係数の**安定性**(`k`
が十分大きいとき、与えられた次数の係数がそれ以上変化しない)を確立し、
極限 `LubinTateEndo` を係数関数として直接定義し、関数等式
`f.subst(LubinTateEndo) = LubinTateEndo.subst(g)` を **exact に**満たす
ことを示す——2変数版 `Found/PGC/LubinTateLimit.lean::LubinTateF` の
1変数版。これで節目(1)(`ResearchPaper/blocked-leaves.json` の
`progress2026_09_04q`〜`s`)の**1変数存在定理そのもの**が完成する:
`g:=f`・`a:=π` の場合には `hf1` から `f` 自身が唯一の解になる
(1変数一意性補題 `powerSeries_uniqueness` と合わせれば分かる)ので、
`[π]_f = f` という既知の事実と矛盾しない。

## 2変数版との違い

`PowerSeries.order` と `MvPowerSeries.order` は defeq でない(§9-1342 で
発見・訂正した点)ため、随所で `PowerSeries.order_eq_order`(mathlib既存の
橋渡し補題)を挟む必要があった——数学的な内容は2変数版と完全に並行。
-/

namespace ABC3.Found.PGC

/-- `φSeq_endo` の逐次関係を露出させる: `φSeq_endo (k+1)` は `φSeq_endo k`
に、`c•X^{k+2}`(`c` はスカラー)を足したもの。 -/
theorem φSeq_endo_succ_eq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (k : ℕ) :
    ∃ c : A,
      (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a (k + 1)).1 =
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 + c • PowerSeries.X ^ (k + 2) := by
  set prev := φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k with hprev_def
  set cex := obstructionVanishesUpTo_step_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf
    prev.2.2 (Nat.succ_ne_zero k) prev.2.1 with hcex_def
  refine ⟨cex.choose, ?_⟩
  show (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a (k + 1)).1 =
    prev.1 + cex.choose • PowerSeries.X ^ (k + 1 + 1)
  rfl

/-- `φSeq_endo` の差の次数は、離れているほど高い: `k≤m` のとき
`(φSeq_endo m) − (φSeq_endo k)` の次数は `≥k+2`。 -/
theorem order_diff_φSeq_endo_ge {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A)
    (k : ℕ) : ∀ m, k ≤ m →
    ((k + 2 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order ((φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 : PowerSeries A) := by
  intro m hkm
  induction m, hkm using Nat.le_induction with
  | base => simp [MvPowerSeries.order_zero]
  | succ n hn ih =>
    obtain ⟨c, hceq⟩ := φSeq_endo_succ_eq hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n
    have hstep : (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a (n + 1)).1 -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 =
        ((φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n).1 -
          (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1) +
          c • (PowerSeries.X : PowerSeries A) ^ (n + 2) := by
      rw [hceq]; ring
    rw [hstep]
    have hcorder : ((n + 2 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 2)) := by
      have hXord : MvPowerSeries.order ((PowerSeries.X : PowerSeries A) ^ (n + 2)) = ((n + 2 : ℕ) : ℕ∞) := by
        rw [← PowerSeries.order_eq_order]
        exact PowerSeries.order_eq.mpr ⟨fun i hi => by
            have hcx := PowerSeries.coeff_X_pow (R := A) i (n + 2)
            rw [if_pos (by exact_mod_cast hi)] at hcx
            rw [hcx]; exact one_ne_zero, fun i hi => by
            rw [PowerSeries.coeff_X_pow]
            rw [if_neg (by intro heq; rw [heq] at hi; exact absurd hi (lt_irrefl _))]⟩
      calc ((n + 2 : ℕ) : ℕ∞) = MvPowerSeries.order ((PowerSeries.X : PowerSeries A) ^ (n + 2)) := hXord.symm
        _ ≤ _ := MvPowerSeries.le_order_smul
    have hnk : ((k + 2 : ℕ) : ℕ∞) ≤ ((n + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : k + 2 ≤ n + 2)
    have hmin : ((k + 2 : ℕ) : ℕ∞) ≤
        min (MvPowerSeries.order ((φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n).1 -
              (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 : PowerSeries A))
          (MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 2))) := by
      rw [le_min_iff]; exact ⟨ih, le_trans hnk hcorder⟩
    exact le_trans hmin MvPowerSeries.min_order_le_add

/-- **係数の安定性**: `e ≤ k+1` のとき、`φSeq_endo m`(`m≥k`)の `e` 係数は
`φSeq_endo k` のそれと一致する。 -/
theorem coeff_φSeq_endo_stable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A)
    (e : ℕ) (k m : ℕ) (hkm : k ≤ m) (he : e ≤ k + 1) :
    PowerSeries.coeff e (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 =
      PowerSeries.coeff e (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 := by
  have hord := order_diff_φSeq_endo_ge hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k m hkm
  have hlt : ((e : ℕ) : ℕ∞) <
      MvPowerSeries.order ((φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 : PowerSeries A) := by
    calc ((e : ℕ) : ℕ∞) ≤ ((k + 1 : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((k + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : k + 1 < k + 2)
      _ ≤ _ := hord
  have hlt' : ((e : ℕ) : ℕ∞) <
      PowerSeries.order ((φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a k).1 : PowerSeries A) := by
    rw [PowerSeries.order_eq_order]; exact hlt
  have hz := PowerSeries.coeff_of_lt_order e hlt'
  rw [map_sub] at hz
  exact sub_eq_zero.mp hz

/-- ★★★**`𝒪_K` 作用の拡張(1変数版)**。`φSeq_endo` の極限——次数 `n` の
係数は、`n` 次まで進んだ近似 `φSeq_endo n` から直接読み取る
(`coeff_φSeq_endo_stable` により、それより先の近似でも同じ値になる)。 -/
noncomputable def LubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries A :=
  PowerSeries.mk (fun n => PowerSeries.coeff n (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n).1)

theorem coeff_LubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (n : ℕ) :
    PowerSeries.coeff n (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) =
      PowerSeries.coeff n (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n).1 :=
  PowerSeries.coeff_mk n _

theorem coeff_LubinTateEndo_eq_φSeq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A)
    (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) =
      PowerSeries.coeff n (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 := by
  rw [coeff_LubinTateEndo]
  exact (coeff_φSeq_endo_stable hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n n m hn (by omega)).symm

theorem order_diff_LubinTateEndo_φSeq_ge {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (m : ℕ) :
    ((m + 1 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 : PowerSeries A) := by
  rw [← PowerSeries.order_eq_order]
  apply PowerSeries.nat_le_order
  intro n hn
  have hnm : n ≤ m := by omega
  have heq := coeff_LubinTateEndo_eq_φSeq hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n m hnm
  rw [map_sub, heq, sub_self]

/-- **g側の安定性**: `LubinTateEndo` を `g` で代入した結果は、`n` 以上まで
進んだ `φSeq_endo m` を代入した結果と、次数 `n` で一致する。 -/
theorem coeff_gsubst_LubinTateEndo_eq_φSeq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A)
    (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n (PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a)) =
      PowerSeries.coeff n (PowerSeries.subst g (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1) := by
  have hHSg : PowerSeries.HasSubst g := by
    show IsNilpotent (PowerSeries.constantCoeff g); rw [hg0]; exact IsNilpotent.zero
  have hsub : PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) -
      PowerSeries.subst g (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 =
      PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1) :=
    (PowerSeries.subst_sub hHSg _ _).symm
  have hdiffge := order_diff_LubinTateEndo_φSeq_ge hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m
  have hdiffge' : ((m + 1 : ℕ) : ℕ∞) ≤
      PowerSeries.order (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 : PowerSeries A) := by
    rw [PowerSeries.order_eq_order]; exact hdiffge
  have horder : ((m + 1 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order (PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1) : PowerSeries A) := by
    have hle := PowerSeries.le_order_subst (g : MvPowerSeries Unit A) hHSg
      (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1)
    have hgorder : (1 : ℕ∞) ≤ MvPowerSeries.order (g : PowerSeries A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hg0
    calc ((m + 1 : ℕ) : ℕ∞) = 1 * ((m + 1 : ℕ) : ℕ∞) := by ring
      _ ≤ MvPowerSeries.order (g : PowerSeries A) * ((m + 1 : ℕ) : ℕ∞) := by gcongr
      _ ≤ MvPowerSeries.order (g : PowerSeries A) * PowerSeries.order
            (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
              (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 : PowerSeries A) := by gcongr
      _ ≤ _ := hle
  have hz : PowerSeries.coeff n
      (PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a -
        (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1)) = 0 := by
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((n : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast hn
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := horder
  rw [← hsub, map_sub] at hz
  exact sub_eq_zero.mp hz

/-- **f側の安定性**: `f` に `LubinTateEndo` を代入した結果は、`n` 以上まで
進んだ `φSeq_endo m` を代入した結果と、次数 `n` で一致する。 -/
theorem coeff_fsubst_LubinTateEndo_eq_φSeq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A)
    (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n (PowerSeries.subst (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) f) =
      PowerSeries.coeff n (PowerSeries.subst (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 f) := by
  set φm := (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).1 with hφm_def
  set δ := LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a - φm with hδ_def
  have hφmadd : φm + δ = LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a := by rw [hδ_def]; ring
  have hφm0 : PowerSeries.constantCoeff φm = 0 := (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).2.2
  have hφmorder : (1 : ℕ∞) ≤ MvPowerSeries.order (φm : PowerSeries A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφm0
  have hδorder : ((m + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (δ : PowerSeries A) :=
    order_diff_LubinTateEndo_φSeq_ge hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m
  have hδ0 : PowerSeries.constantCoeff δ = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((0 : ℕ) : ℕ∞) < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : 0 < m + 1)
      _ ≤ _ := hδorder
  have hlin := coeff_subst_linearize_1var hφm0 hδ0 hφmorder hδorder (by omega : 1 ≤ m + 1) f π hf0 hf1 n
    (by exact_mod_cast (by omega : n ≤ m + 1))
  rw [hφmadd] at hlin
  have hδn0 : PowerSeries.coeff n δ = 0 := by
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((n : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast hn
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := hδorder
  rw [hδn0, mul_zero] at hlin
  exact sub_eq_zero.mp hlin

theorem coeff_LubinTateEndo_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (n : ℕ) :
    PowerSeries.coeff n (PowerSeries.subst (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) f) =
      PowerSeries.coeff n (PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a)) := by
  set m := n with hm_def
  have h1 := coeff_fsubst_LubinTateEndo_eq_φSeq hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n m le_rfl
  have h2 := coeff_gsubst_LubinTateEndo_eq_φSeq hq hπmax g hg0 hg1 hg f hf0 hf1 hf a n m le_rfl
  have h3 := (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a m).2.1 n (by omega)
  rw [map_sub] at h3
  rw [h1, h2]
  exact sub_eq_zero.mp h3

/-- ★★★★★★★★**`𝒪_K` 作用の拡張(1変数版)の存在補題**。`LubinTateEndo` は
関数等式 `f.subst(LubinTateEndo) = LubinTateEndo.subst(g)`(すなわち
`f(φ(X))=φ(g(X))`)を power series の等式として満たす——2変数版
`LubinTateF_functional_equation` の1変数版。 -/
theorem LubinTateEndo_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.subst (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) f =
      PowerSeries.subst g (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) :=
  PowerSeries.ext (coeff_LubinTateEndo_functional_equation hq hπmax g hg0 hg1 hg f hf0 hf1 hf a)

theorem constantCoeff_LubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.constantCoeff (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_LubinTateEndo,
    PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a 0).2.2

/-- `LubinTateEndo` の次数1の係数は目標どおり `a`——`φSeq_endo` の出発点
`a•X` から、次数2以上の補正しか足されないから。 -/
theorem coeff_one_LubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.coeff 1 (LubinTateEndo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a) = a := by
  rw [coeff_LubinTateEndo]
  obtain ⟨c, hceq⟩ := φSeq_endo_succ_eq hq hπmax g hg0 hg1 hg f hf0 hf1 hf a 0
  show PowerSeries.coeff 1 (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a 1).1 = a
  rw [hceq]
  have hbase : (φSeq_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf a 0).1 = a • (PowerSeries.X : PowerSeries A) := rfl
  rw [hbase, map_add, PowerSeries.coeff_smul, PowerSeries.coeff_X, if_pos rfl, smul_eq_mul, mul_one]
  rw [PowerSeries.coeff_smul, PowerSeries.coeff_X_pow, if_neg (by omega), smul_zero, add_zero]

end ABC3.Found.PGC
