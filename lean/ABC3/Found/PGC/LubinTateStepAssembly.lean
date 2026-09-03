import ABC3.Found.PGC.LubinTateGLinearization
import ABC3.Found.PGC.LubinTateStep

/-!
# Lubin-Tate 形式群の存在補題: 次数ごとの帰納法1ステップの組み立て(`sorry` 無し)

これまでに確立した3部品——可除性(`Found/PGC/LubinTateStep.lean::
exists_step_solution_for_R`)、f側の線形化(`Found/PGC/LubinTateLinearization.
lean::coeff_subst_linearize`)、g側の線形化(`Found/PGC/LubinTateGLinearization.
lean::coeff_subst_g_linearize`)——を実際に組み合わせて、「`Φ` の障害
`Obstruction(Φ) := Φ.subst(g,g) − f.subst(Φ)` が次数 `n` まで消えている」
→「`Φ' := Φ+φ`(`φ` は次数 `n+1` の斉次式)の障害が次数 `n+1` まで消えて
いる」という**帰納法の1ステップ全体**を確立する。

## 鍵になった簡略化: `φ` の斉次性は「取り出す」だけで済む

`exists_step_solution_for_R` が返す `φ₀` はそれ自体が次数 `n+1` の斉次式で
あるとは限らない(係数ごとの独立な選択、`Classical.choice` 経由)。しかし
`coeff_subst_g_linearize` は `φ` が**厳密に**斉次であることを要求する
(f側の `coeff_subst_linearize` は `φ.order ≥ n+1` だけで足りるが、g側は
それでは足りない)。

ここで `φ := MvPowerSeries.homogeneousComponent (n+1) φ₀` と**取り出す**だけで
両方の要求が満たせることに気づいた: `homogeneousComponent (n+1)` は
`R`-線型写像(`→ₗ[R]`)なので `(π−π^{n+1})•φ₀ = R_n` の両辺に適用すると
`(π−π^{n+1})•φ = homogeneousComponent(n+1)(R_n) = R_n`(`R_n` 自身が既に
次数 `n+1` の斉次式なので、`homogeneousComponent(n+1)` は不動点——
`isHomogeneous_iff_eq_homogeneousComponent` の言い換え)。しかも
`φ` の斉次性そのものは `isHomogeneous_homogeneousComponent` から無料で
出る。「整域でのキャンセル」のような別の議論は一切不要だった。

## まだ無いもの

本ファイルは**1ステップ**の帰納法補題を確立する。実際に `Nat.rec` で
`Φ : ℕ → MvPowerSeries (Fin 2) A` の無限列を構成し、極限 `F` を
`MvPowerSeries.mk`(または `homogeneousComponent` の無限和の代替)で
組み立て、`F` が関数等式 `F.subst(g,g) = f.subst(F)` を**exact に**満たす
ことを示す段は、まだ残る——本ファイルの1ステップ補題を `n` について
`Nat.rec` で繰り返し適用するだけで進むが、「無限次までの整合性」から
「実際の等式」への橋渡し(`MvPowerSeries` は係数関数そのものなので、
各次数の係数が最終的に安定することを言えば `F` の各係数が定まる)は
別途の組み立て作業として残る。
-/

namespace ABC3.Found.PGC

open MvPowerSeries in
/-- ★★★**次数ごとの帰納法の1ステップ**。`Obstruction Φ := Φ.subst(g,g) −
f.subst(Φ)` が次数 `≤n` の範囲で消えているとき、次数 `n+1` の斉次式 `φ` が
存在して、`Obstruction (Φ+φ)` は次数 `≤n+1` の範囲で消える。3部品
(可除性・f側線形化・g側線形化)の組み合わせで、φ を `exists_step_solution_
for_R` の解の次数 `n+1` 成分として取ることで得られる。 -/
theorem exists_next_step {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    {Φ : MvPowerSeries (Fin 2) A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    {n : ℕ} (hn : n ≠ 0)
    (hinv : ∀ e : Fin 2 →₀ ℕ, Finsupp.degree e ≤ n →
      MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
          PowerSeries.subst Φ f) = 0) :
    ∃ φ : MvPowerSeries (Fin 2) A,
      (∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ = 0) ∧
      ∀ e : Fin 2 →₀ ℕ, Finsupp.degree e ≤ n + 1 →
        MvPowerSeries.coeff e
          (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) (Φ + φ) -
            PowerSeries.subst (Φ + φ) f) = 0 := by
  obtain ⟨φ₀, hφ₀eq⟩ := exists_step_solution_for_R hq hπmax g hg0 f hf hg Φ hΦ0 n hn
  set φ : MvPowerSeries (Fin 2) A := MvPowerSeries.homogeneousComponent (n + 1) φ₀ with hφ_def
  have hφHomog : MvPowerSeries.IsHomogeneous φ (n + 1) :=
    MvPowerSeries.isHomogeneous_homogeneousComponent φ₀ (n + 1)
  have hφ : ∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ = 0 :=
    fun d hd => hφHomog.coeff_eq_zero hd
  set Rn : MvPowerSeries (Fin 2) A :=
    MvPowerSeries.homogeneousComponent (n + 1)
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
        PowerSeries.subst Φ f) with hRn_def
  have hRnHomog : MvPowerSeries.IsHomogeneous Rn (n + 1) :=
    MvPowerSeries.isHomogeneous_homogeneousComponent _ (n + 1)
  have hRnFix : MvPowerSeries.homogeneousComponent (n + 1) Rn = Rn :=
    (MvPowerSeries.isHomogeneous_iff_eq_homogeneousComponent.mp hRnHomog).symm
  have hφeq : (π - π ^ (n + 1)) • φ = Rn := by
    rw [hφ_def]
    calc (π - π ^ (n + 1)) • MvPowerSeries.homogeneousComponent (n + 1) φ₀
        = MvPowerSeries.homogeneousComponent (n + 1) ((π - π ^ (n + 1)) • φ₀) :=
          (map_smul (MvPowerSeries.homogeneousComponent (n + 1)) (π - π ^ (n + 1)) φ₀).symm
      _ = MvPowerSeries.homogeneousComponent (n + 1) Rn := by rw [hφ₀eq]
      _ = Rn := hRnFix
  have hφ0 : MvPowerSeries.constantCoeff φ = 0 := by
    rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact hφ (0 : Fin 2 →₀ ℕ) (by simp)
  have hφord : ((n + 1 : ℕ) : ℕ∞) ≤ φ.order := by
    apply MvPowerSeries.nat_le_order
    intro d hd
    apply hφ
    intro hcontra
    rw [hcontra] at hd
    exact absurd hd (lt_irrefl _)
  have hΦorder : 1 ≤ Φ.order := MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hΦ0
  have haHS := hasSubst_g_subst_X (σ := Fin 2) g hg0
  refine ⟨φ, hφ, fun e he => ?_⟩
  have hAdd : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) (Φ + φ) =
      MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ +
        MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) φ :=
    MvPowerSeries.subst_add haHS Φ φ
  have hGside : MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) φ) =
      π ^ (n + 1) * MvPowerSeries.coeff e φ :=
    coeff_subst_g_linearize hφ π hg0 hg1 e he
  have hFside : MvPowerSeries.coeff e (PowerSeries.subst (Φ + φ) f) -
      MvPowerSeries.coeff e (PowerSeries.subst Φ f) = π * MvPowerSeries.coeff e φ :=
    coeff_subst_linearize hΦ0 hφ0 hΦorder hφord (by omega) f π hf0 hf1 e (by exact_mod_cast he)
  have hkey : MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) (Φ + φ) -
        PowerSeries.subst (Φ + φ) f) =
      MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
          PowerSeries.subst Φ f) - (π - π ^ (n + 1)) * MvPowerSeries.coeff e φ := by
    rw [map_sub, hAdd, map_add, map_sub]
    linear_combination hGside - hFside
  rw [hkey]
  rcases (by omega : Finsupp.degree e ≤ n ∨ Finsupp.degree e = n + 1) with hle | heq
  · rw [hinv e hle, hφ e (by omega)]; ring
  · have hcoeffφeq : (π - π ^ (n + 1)) * MvPowerSeries.coeff e φ =
        MvPowerSeries.coeff e
          (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
            PowerSeries.subst Φ f) := by
      have hcoeffRn : MvPowerSeries.coeff e Rn =
          MvPowerSeries.coeff e
            (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
              PowerSeries.subst Φ f) := by
        rw [hRn_def, MvPowerSeries.coeff_homogeneousComponent, if_pos heq]
      rw [← hcoeffRn, ← hφeq, MvPowerSeries.coeff_smul]
    rw [hcoeffφeq]; ring

end ABC3.Found.PGC
