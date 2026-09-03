import ABC3.Found.PGC.LubinTateBaseCase

/-!
# Lubin-Tate 形式群の存在補題: 次数ごとの近似列(`sorry` 無し)

`Found/PGC/LubinTateBaseCase.lean`(出発点、次数1まで解けている `X_0+X_1`)と
`Found/PGC/LubinTateStepAssembly.lean::exists_next_step`(1ステップ)を
`Nat.rec` で繋いで、**任意の次数まで解けている近似**の列
`ΦSeq k : MvPowerSeries (Fin 2) A`(`k`次まで `Obstruction` が消えている)を
構成する。

## まだ無いもの

この列の**極限** `F`(各係数が `ΦSeq k` の対応する係数で `k` が十分大きい
とき安定する、という事実を使って係数関数として直接定義する)を構成し、
`F` が関数等式 `F.subst(g,g) = f.subst(F)` を**exact に**満たすことを示す
段は、まだ残る——`MvPowerSeries.truncTotal_subst_*` 系(代入が引数の低次の
切り捨てにしか依存しないという「連続性」)がこの段の土台になりうると
見ている。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hgres : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

/-- `Φ` の障害(`Φ.subst(g,g) − f.subst(Φ)`)が次数 `≤n` の範囲で消えている、
という次数ごとの再帰の不変量。 -/
def ObstructionVanishesUpTo (Φ : MvPowerSeries (Fin 2) A) (n : ℕ) : Prop :=
  ∀ e : Fin 2 →₀ ℕ, Finsupp.degree e ≤ n →
    MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
        PowerSeries.subst Φ f) = 0

include hg0 hg1 hf0 hf1 in
omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- 出発点 `X_0+X_1` は次数1まで解けている。 -/
theorem obstructionVanishesUpTo_base :
    ObstructionVanishesUpTo g f ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) 1 :=
  fun e he => base_case_invariant g hg0 hg1 f hf0 hf1 e he

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- 1ステップ: 次数 `n`(`n≠0`)まで解けている `Φ` から、次数 `n+1` まで
解けている `Φ+φ`(`φ` は次数 `n+1` の斉次式)が作れる。 -/
theorem obstructionVanishesUpTo_step {Φ : MvPowerSeries (Fin 2) A}
    (hΦ0 : MvPowerSeries.constantCoeff Φ = 0) {n : ℕ} (hn : n ≠ 0)
    (hinv : ObstructionVanishesUpTo g f Φ n) :
    ∃ φ : MvPowerSeries (Fin 2) A,
      (∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ = 0) ∧
      ObstructionVanishesUpTo g f (Φ + φ) (n + 1) :=
  exists_next_step hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres hΦ0 hn hinv

/-- 次数ごとの近似列が満たすべき不変量: 障害が次数 `k+1` まで消えていて、
かつ定数項が0(次のステップを踏むために要る)。 -/
def SeqInvariant (Φ : MvPowerSeries (Fin 2) A) (k : ℕ) : Prop :=
  ObstructionVanishesUpTo g f Φ (k + 1) ∧ MvPowerSeries.constantCoeff Φ = 0

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- ★★**次数ごとの近似列**。`ΦSeq k` は次数 `k+1` まで `Obstruction` が
消えている `MvPowerSeries (Fin 2) A`(定数項0)。`Nat.rec` で `Found/PGC/
LubinTateBaseCase.lean` の出発点から `Found/PGC/LubinTateStepAssembly.lean::
exists_next_step` を繰り返し適用して作る。 -/
noncomputable def ΦSeq : (k : ℕ) → {Φ : MvPowerSeries (Fin 2) A // SeqInvariant g f Φ k} :=
  fun k => Nat.rec
    (⟨(MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1,
        obstructionVanishesUpTo_base g hg0 hg1 f hf0 hf1, by
          simp [map_add, MvPowerSeries.constantCoeff_X]⟩ :
      {Φ : MvPowerSeries (Fin 2) A // SeqInvariant g f Φ 0})
    (fun k prev =>
      let φex := obstructionVanishesUpTo_step hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres
        prev.2.2 (Nat.succ_ne_zero k) prev.2.1
      have hφ0 : MvPowerSeries.constantCoeff φex.choose = 0 := by
        rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
        exact φex.choose_spec.1 (0 : Fin 2 →₀ ℕ) (by simp)
      (⟨prev.1 + φex.choose, φex.choose_spec.2, by
          rw [map_add, prev.2.2, hφ0, add_zero]⟩ :
        {Φ : MvPowerSeries (Fin 2) A // SeqInvariant g f Φ (k + 1)}))
    k

end ABC3.Found.PGC
