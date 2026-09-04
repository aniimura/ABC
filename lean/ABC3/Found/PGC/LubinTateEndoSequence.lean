import ABC3.Found.PGC.LubinTateEndoBaseCase

/-!
# `𝒪_K` 作用への拡張(1変数版)の次数ごとの近似列(`sorry` 無し)

`Found/PGC/LubinTateEndoBaseCase.lean`(出発点、次数1まで解けている `a•X`)と
`Found/PGC/LubinTateEndoStepAssembly.lean::exists_next_step_endo`(1ステップ)
を `Nat.rec` で繋いで、**任意の次数まで解けている近似**の列
`φSeq_endo k : PowerSeries A`(`k+1`次まで `Obstruction` が消えている)を
構成する——2変数版 `LubinTateSequence.lean::ΦSeq` の1変数版。

## まだ無いもの

この列の**極限** `φ_∞`(各係数が `φSeq_endo k` の対応する係数で `k` が
十分大きいとき安定する、という事実を使って係数関数として直接定義する)を
構成し、`φ_∞` が関数等式 `f.subst(φ_∞) = φ_∞.subst(g)` を**exact に**
満たすことを示す段は、2変数版(`LubinTateLimit.lean`)と同型だが1変数へ
作り直す必要があり、別途の課題として残る。
-/

namespace ABC3.Found.PGC

/-- `φ` の障害(`f.subst(φ) − g.subst(φ)`)が次数 `≤n` の範囲で消えている、
という次数ごとの再帰の不変量。 -/
def ObstructionVanishesUpTo_endo {A : Type*} [CommRing A] (g f : PowerSeries A)
    (φ : PowerSeries A) (n : ℕ) : Prop :=
  ∀ k ≤ n, PowerSeries.coeff k (PowerSeries.subst φ f - PowerSeries.subst g φ) = 0

/-- 出発点 `a•X` は次数1まで解けている。 -/
theorem obstructionVanishesUpTo_base_endo {A : Type*} [CommRing A] {π a : A}
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π) :
    ObstructionVanishesUpTo_endo g f (a • (PowerSeries.X : PowerSeries A)) 1 :=
  fun k hk => base_case_invariant_endo g hg0 hg1 f hf0 hf1 k hk

/-- 1ステップ: 次数 `n`(`n≠0`)まで解けている `φ` から、次数 `n+1` まで
解けている `φ+c•X^{n+1}`(`c` はスカラー)が作れる。 -/
theorem obstructionVanishesUpTo_step_endo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    {φ : PowerSeries A} (hφ0 : PowerSeries.constantCoeff φ = 0)
    {n : ℕ} (hn : n ≠ 0)
    (hinv : ObstructionVanishesUpTo_endo g f φ n) :
    ∃ c : A, PowerSeries.constantCoeff (φ + c • PowerSeries.X ^ (n + 1)) = 0 ∧
      ObstructionVanishesUpTo_endo g f (φ + c • PowerSeries.X ^ (n + 1)) (n + 1) :=
  exists_next_step_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf hφ0 hn hinv

/-- 次数ごとの近似列が満たすべき不変量: 障害が次数 `k+1` まで消えていて、
かつ定数項が0(次のステップを踏むために要る)。 -/
def SeqInvariant_endo {A : Type*} [CommRing A] (g f : PowerSeries A) (φ : PowerSeries A)
    (k : ℕ) : Prop :=
  ObstructionVanishesUpTo_endo g f φ (k + 1) ∧ PowerSeries.constantCoeff φ = 0

/-- ★★**次数ごとの近似列(1変数版)**。`φSeq_endo k` は次数 `k+1` まで
`Obstruction` が消えている `PowerSeries A`(定数項0)。`Nat.rec` で
`Found/PGC/LubinTateEndoBaseCase.lean` の出発点から `Found/PGC/
LubinTateEndoStepAssembly.lean::exists_next_step_endo` を繰り返し適用して
作る——2変数版 `ΦSeq` と全く同じ構造だが、各段で足す「補正」がスカラー
`c` そのもの(`homogeneousComponent` の抽出が不要)なぶん単純になった。 -/
noncomputable def φSeq_endo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    (k : ℕ) → {φ : PowerSeries A // SeqInvariant_endo g f φ k} :=
  fun k => Nat.rec
    (⟨a • (PowerSeries.X : PowerSeries A),
        obstructionVanishesUpTo_base_endo g hg0 hg1 f hf0 hf1, by
          show PowerSeries.constantCoeff (a • (PowerSeries.X : PowerSeries A)) = 0
          show MvPowerSeries.constantCoeff (a • (PowerSeries.X : PowerSeries A)) = 0
          rw [MvPowerSeries.constantCoeff_smul]
          show a • MvPowerSeries.constantCoeff (MvPowerSeries.X () : MvPowerSeries Unit A) = 0
          rw [MvPowerSeries.constantCoeff_X, smul_zero]⟩ :
      {φ : PowerSeries A // SeqInvariant_endo g f φ 0})
    (fun k prev =>
      let cex := obstructionVanishesUpTo_step_endo hq hπmax g hg0 hg1 hg f hf0 hf1 hf
        prev.2.2 (Nat.succ_ne_zero k) prev.2.1
      (⟨prev.1 + cex.choose • (PowerSeries.X : PowerSeries A) ^ (k + 1 + 1),
          cex.choose_spec.2, cex.choose_spec.1⟩ :
        {φ : PowerSeries A // SeqInvariant_endo g f φ (k + 1)}))
    k

end ABC3.Found.PGC
