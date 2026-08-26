/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithFunctor

/-!
# 算術因子の引き戻しの関手性 —— 塔での推移(`Example 6.3` の `pull_comp` / `pull_id`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★`ArithDatum` が要求する 2 本

`ArithFrobenioid.lean` の `ArithDatum` は `pull_id` と `pull_comp` を要求する。
★`arithExtend`(`ArithFunctor.lean`)がそれを満たすことを、ここで示す。

## ★★★中身は「局所次数は塔で掛け算になる」

  `[N_V : L_v] = [M_W : L_v] · [N_V : M_W]`   (`W = V|_M`, `v = V|_L`)

| 成分 | 根拠 |
|---|---|
| 非アルキメデス `e` | `Ideal.ramificationIdx_algebra_tower'` |
| 非アルキメデス `f` | `Ideal.inertiaDeg_algebra_tower` |
| アルキメデス `mult` | `mult` は `1` か `2` で、塔に沿って割り切る |

★素点の制限が推移すること(`resPlace_resPlace`)は
`Ideal.under_under`(非アルキメデス)と `IsScalarTower.algebraMap_eq`(アルキメデス)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `resHOS_resHOS` / `resFin_resFin` / `resInf_resInf` / `resPlace_resPlace` | 素点の制限は推移する |
| `ramIdx_tower` / `inertDeg_tower` / `localDeg_tower` | ★局所次数は塔で掛け算 |
| `arithExtend_comp` | ★★★**`pull_comp`** |
| `arithExtend_id` | ★**`pull_id`** |
-/

namespace ABC3.Found.Divisor

open NumberField IsDedekindDomain Ideal

universe u

section Tower

variable {L M N : Type u} [Field L] [Field M] [Field N]
  [NumberField L] [NumberField M] [NumberField N]
  [Algebra L M] [Algebra M N] [Algebra L N] [IsScalarTower L M N]

/-! ## ★1. 素点の制限は推移する -/

omit [NumberField L] [NumberField M] [NumberField N] in
/-- ★高さ 1 の素点の制限は推移する(`Ideal.under_under`)。 -/
theorem resHOS_resHOS (V : HeightOneSpectrum (𝓞 N)) :
    resHOS (L := L) (resHOS (L := M) V) = resHOS (L := L) V := by
  refine HeightOneSpectrum.ext ?_
  exact Ideal.under_under V.asIdeal

/-- ★有限素点の制限は推移する。 -/
theorem resFin_resFin (W : FinitePlace N) :
    resFin (L := L) (resFin (L := M) W) = resFin (L := L) W := by
  rw [resFin, maximalIdeal_resFin, resHOS_resHOS, resFin]

omit [NumberField L] [NumberField M] [NumberField N] in
/-- ★無限素点の制限は推移する。 -/
theorem resInf_resInf (W : InfinitePlace N) :
    resInf (L := L) (resInf (L := M) W) = resInf (L := L) W := by
  rw [resInf, resInf, resInf, IsScalarTower.algebraMap_eq L M N,
    InfinitePlace.comap_comp]

/-- ★★**素点の制限は推移する**。 -/
theorem resPlace_resPlace (V : ArithPlace N) :
    resPlace (L := L) (resPlace (L := M) V) = resPlace (L := L) V := by
  cases V with
  | inl W => simp [resFin_resFin]
  | inr W => simp [resInf_resInf]

/-! ## ★2. 局所次数は塔で掛け算 -/

omit [NumberField L] in
/-- ★★**分岐指数は塔で掛け算**。 -/
theorem ramIdx_tower (W : FinitePlace N) :
    ramIdx (L := L) W = ramIdx (L := L) (resFin (L := M) W) * ramIdx (L := M) W := by
  haveI : (W.maximalIdeal.asIdeal).LiesOver ((resHOS (L := M) W.maximalIdeal).asIdeal) :=
    liesOver_resHOS W.maximalIdeal
  haveI : ((resHOS (L := M) W.maximalIdeal).asIdeal).LiesOver
      ((resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).asIdeal) :=
    liesOver_resHOS (L := L) (resHOS (L := M) W.maximalIdeal)
  rw [ramIdx, ramIdx, ramIdx, maximalIdeal_resFin, resHOS_resHOS]
  have h := Ideal.ramificationIdx_algebra_tower' (R := 𝓞 L) (S := 𝓞 M) (T := 𝓞 N)
    (resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).asIdeal
    (resHOS (L := M) W.maximalIdeal).asIdeal W.maximalIdeal.asIdeal
  rw [resHOS_resHOS] at h
  exact h

/-- ★★**剰余次数は塔で掛け算**。 -/
theorem inertDeg_tower (W : FinitePlace N) :
    inertDeg (L := L) W = inertDeg (L := L) (resFin (L := M) W) * inertDeg (L := M) W := by
  haveI : (W.maximalIdeal.asIdeal).LiesOver ((resHOS (L := M) W.maximalIdeal).asIdeal) :=
    liesOver_resHOS W.maximalIdeal
  haveI : ((resHOS (L := M) W.maximalIdeal).asIdeal).LiesOver
      ((resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).asIdeal) :=
    liesOver_resHOS (L := L) (resHOS (L := M) W.maximalIdeal)
  haveI hmL : ((resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).asIdeal).IsMaximal :=
    Ideal.IsPrime.isMaximal (resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).isPrime
      (resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).ne_bot
  haveI hmM : ((resHOS (L := M) W.maximalIdeal).asIdeal).IsMaximal :=
    Ideal.IsPrime.isMaximal (resHOS (L := M) W.maximalIdeal).isPrime
      (resHOS (L := M) W.maximalIdeal).ne_bot
  rw [inertDeg, inertDeg, inertDeg, maximalIdeal_resFin, resHOS_resHOS]
  have h := Ideal.inertiaDeg_algebra_tower (R := 𝓞 L) (S := 𝓞 M) (T := 𝓞 N)
    (resHOS (L := L) (resHOS (L := M) W.maximalIdeal)).asIdeal
    (resHOS (L := M) W.maximalIdeal).asIdeal W.maximalIdeal.asIdeal
  rw [resHOS_resHOS] at h
  exact h

/-- ★★★**局所次数は塔で掛け算**。

★非アルキメデスは `e·f` の掛け算、アルキメデスは `mult` が `1` か `2` で
塔に沿って割り切ることから。 -/
theorem localDeg_tower (V : ArithPlace N) :
    localDeg (L := L) V = localDeg (L := L) (resPlace (L := M) V) * localDeg (L := M) V := by
  cases V with
  | inl W =>
      rw [localDeg_inl, resPlace_inl, localDeg_inl, localDeg_inl,
        ramIdx_tower (L := L) (M := M) W, inertDeg_tower (L := L) (M := M) W]
      ring
  | inr W =>
      rw [localDeg_inr, resPlace_inr, localDeg_inr, localDeg_inr, resInf_resInf]
      have h1 : (resInf (L := L) W).mult ∣ (resInf (L := M) W).mult := by
        have := mult_resInf_dvd (L := L) (M := M) (resInf (L := M) W)
        rwa [resInf_resInf] at this
      have h2 : (resInf (L := M) W).mult ∣ W.mult := mult_resInf_dvd (L := M) W
      have h3 : (resInf (L := L) W).mult ∣ W.mult := mult_resInf_dvd (L := L) W
      rcases mult_eq_one_or_two (resInf (L := L) W) with e1 | e1 <;>
        rcases mult_eq_one_or_two (resInf (L := M) W) with e2 | e2 <;>
        rcases mult_eq_one_or_two W with e3 | e3 <;>
        simp only [e1, e2, e3] at h1 h2 h3 ⊢ <;> omega

/-! ## ★3. `pull_comp` -/

/-- ★★★★★**引き戻しは合成と両立する**(`ArithDatum.pull_comp`)。 -/
theorem arithExtend_comp (d : ArithPlace L →₀ ℝ) :
    arithExtend (L := L) (M := N) d
      = arithExtend (L := M) (M := N) (arithExtend (L := L) (M := M) d) := by
  refine Finsupp.ext fun V => ?_
  show (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V)
    = (localDeg (L := M) V : ℝ)
        * (arithExtend (L := L) (M := M) d) (resPlace (L := M) V)
  show (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V)
    = (localDeg (L := M) V : ℝ)
        * ((localDeg (L := L) (resPlace (L := M) V) : ℝ)
            * d (resPlace (L := L) (resPlace (L := M) V)))
  rw [resPlace_resPlace, localDeg_tower (L := L) (M := M) V]
  push_cast
  ring

end Tower

/-! ## ★4. `pull_id` -/

section Id

variable {L : Type u} [Field L] [NumberField L]

omit [NumberField L] in
theorem resHOS_self (V : HeightOneSpectrum (𝓞 L)) : resHOS (L := L) V = V := by
  refine HeightOneSpectrum.ext ?_
  show V.asIdeal.comap (algebraMap (𝓞 L) (𝓞 L)) = V.asIdeal
  simp

theorem resFin_self (W : FinitePlace L) : resFin (L := L) W = W := by
  rw [resFin, resHOS_self, FinitePlace.mk_maximalIdeal]

omit [NumberField L] in
theorem resInf_self (W : InfinitePlace L) : resInf (L := L) W = W := by
  show W.comap (algebraMap L L) = W
  rw [show (algebraMap L L) = RingHom.id L from rfl]
  exact InfinitePlace.comap_id W

theorem resPlace_self (V : ArithPlace L) : resPlace (L := L) V = V := by
  cases V with
  | inl W => simp [resFin_self]
  | inr W => simp [resInf_self]

theorem localDeg_self (V : ArithPlace L) : localDeg (L := L) V = 1 := by
  cases V with
  | inl W =>
      rw [localDeg_inl, ramIdx, inertDeg, resHOS_self]
      haveI : (W.maximalIdeal.asIdeal).IsMaximal :=
        Ideal.IsPrime.isMaximal W.maximalIdeal.isPrime W.maximalIdeal.ne_bot
      haveI : (W.maximalIdeal.asIdeal).LiesOver (W.maximalIdeal.asIdeal) := ⟨rfl⟩
      have hmap : Ideal.map (algebraMap (𝓞 L) (𝓞 L)) W.maximalIdeal.asIdeal
          = W.maximalIdeal.asIdeal := by
        show Ideal.map (RingHom.id (𝓞 L)) W.maximalIdeal.asIdeal = _
        simp
      have hne_top : Ideal.map (algebraMap (𝓞 L) (𝓞 L)) W.maximalIdeal.asIdeal ≠ ⊤ := by
        rw [hmap]; exact W.maximalIdeal.isPrime.ne_top
      have hne_bot : Ideal.map (algebraMap (𝓞 L) (𝓞 L)) W.maximalIdeal.asIdeal ≠ ⊥ := by
        rw [hmap]; exact W.maximalIdeal.ne_bot
      have he := Ideal.ramificationIdx_map_self_eq_one hne_top hne_bot
      rw [hmap] at he
      have hf : Ideal.inertiaDeg W.maximalIdeal.asIdeal W.maximalIdeal.asIdeal = 1 := by
        rw [Ideal.inertiaDeg_algebraMap]
        exact Module.finrank_self _
      show W.maximalIdeal.asIdeal.ramificationIdx W.maximalIdeal.asIdeal
           * W.maximalIdeal.asIdeal.inertiaDeg W.maximalIdeal.asIdeal = 1
      rw [he, hf]
  | inr W =>
      rw [localDeg_inr, resInf_self, Nat.div_self]
      rcases mult_eq_one_or_two W with h | h <;> omega

/-- ★★★★**引き戻しは恒等射で恒等**(`ArithDatum.pull_id`)。 -/
theorem arithExtend_id (d : ArithPlace L →₀ ℝ) : arithExtend (L := L) (M := L) d = d := by
  refine Finsupp.ext fun V => ?_
  show (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V) = d V
  rw [resPlace_self, localDeg_self]
  simp

end Id

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Example 6.3` の因子の引き戻しの関手性(`pull_id` / `pull_comp`)。 -/
def arithExtend_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 因子の引き戻しの関手性(局所次数は塔で掛け算)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
