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

theorem emb_injective : Function.Injective (emb (A := A)) := by
  intro p q hpq
  ext n
  have := coeff_emb (A := A) p (Finsupp.single (0 : Fin 2) n)
  have hq := coeff_emb (A := A) q (Finsupp.single (0 : Fin 2) n)
  simp only [Finsupp.single_eq_same] at this hq
  rw [hpq] at this
  rw [this] at hq
  simpa using hq

end ABC3.Found.PGC
