import ABC3.Found.GaloisRep.CountInfty

/-!
# Galois (G5) 第 188 ブロック —— **★★★★★★★★`∏_{i<n} τ_{iP}^*(f_P)` は定数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★交代性の中核

第 187 で極の位数が決まったので、各素点での指数が**全部書ける**:

    count_v(τ_T f_P) = n·[Q_v + T = P] − n·[Q_v + T = O]

★`T = O` のときは `τ = id` として同じ式が成り立つ(`Q_v ≠ O` だから第 2 項は 0)。

### ★★★★★★和が消える理由——`i ↦ i+1 mod n`

    count_v(∏_{i<n} τ_{iP}^* f_P)
      = Σ_{i<n} n·[Q_v + iP = P] − Σ_{i<n} n·[Q_v + iP = O]

★`i ↦ (i+1) mod n` は `range n` の全単射で、`nP = O` より

    Q_v + ((i+1) mod n)·P = (Q_v + iP) + P

★★したがって「`Q_v + iP = O`」の個数と「`Q_v + iP = P`」の個数は**一致**する。
★★★**`P` の位数が `n` より小さくても構わない**——両辺が同じだけ重複するだけである。

### ★★★★★★指数が全部 0 なら定数

`spanSingleton z = 1` となり、`z` は `F[W]` の単元の像。★第 128 で定数。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `nsmul_mod` | `n•P = 0` なら `k•P` は `k mod n` で決まる |
| `sum_shift_eq` | ★★★★★**`i ↦ i+1 mod n` による個数の一致** |
| `placeOf_eq_iff` | 素点が等しい ⟺ 点が等しい |
| `count_prod'` | `count` は有限積を和にする |
| `const_of_count_eq_zero` | ★★★★★★**指数が全部 0 なら定数** |
| `IsTranslate` | `τ` が点 `T` による平行移動であること(`T = O` 込み) |
| `count_translate_fP_gen` | ★★★★★★★**平行移動した `f_P` の指数(全素点)** |
| `prodTranslate_const` | ★★★★★★★★**積は定数** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

universe u

/-! ## ★数え上げの段 -/

/-- `n • P = 0` なら `k • P` は `k mod n` だけで決まる。 -/
theorem nsmul_mod {A : Type u} [AddCommGroup A] (n : ℕ) {P : A} (hP : n • P = 0) (k : ℕ) :
    (k % n) • P = k • P := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn; simp
  · conv_rhs => rw [← Nat.div_add_mod k n]
    rw [add_smul, mul_comm, mul_smul, hP, smul_zero, zero_add]

theorem succ_mod_roundtrip {n a : ℕ} (hn : 1 ≤ n) (ha : a < n) :
    ((a + 1) % n + (n - 1)) % n = a := by
  rcases Nat.lt_or_ge (a + 1) n with hlt | hge
  · have hm1 : (a + 1) % n = a + 1 := Nat.mod_eq_of_lt hlt
    have he : a + 1 + (n - 1) = a + n := by omega
    rw [hm1, he, Nat.add_mod_right, Nat.mod_eq_of_lt ha]
  · have han : a + 1 = n := by omega
    have hm1 : (a + 1) % n = 0 := by rw [han, Nat.mod_self]
    have he : 0 + (n - 1) = n - 1 := by omega
    rw [hm1, he, Nat.mod_eq_of_lt (by omega : n - 1 < n)]
    omega

theorem pred_mod_roundtrip {n a : ℕ} (hn : 1 ≤ n) (ha : a < n) :
    ((a + (n - 1)) % n + 1) % n = a := by
  rcases Nat.eq_zero_or_pos a with rfl | hpos
  · have he : 0 + (n - 1) = n - 1 := by omega
    have hm1 : (0 + (n - 1)) % n = n - 1 := by rw [he]; exact Nat.mod_eq_of_lt (by omega)
    have he2 : n - 1 + 1 = n := by omega
    rw [hm1, he2, Nat.mod_self]
  · have he : a + (n - 1) = (a - 1) + n := by omega
    have hm1 : (a + (n - 1)) % n = a - 1 := by
      rw [he, Nat.add_mod_right]; exact Nat.mod_eq_of_lt (by omega)
    have he2 : a - 1 + 1 = a := by omega
    rw [hm1, he2, Nat.mod_eq_of_lt ha]

/-- ★★★★★**「`Q + iP = R + P`」と「`Q + iP = R`」の個数は等しい**——`i ↦ i+1 mod n`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`P` の位数が `n` より小さくても、両辺が同じだけ重複するので成り立つ。 -/
theorem sum_shift_eq {A : Type u} [AddCommGroup A] [DecidableEq A] (n : ℕ) (hn : 1 ≤ n)
    {P : A} (hP : n • P = 0) (Q R : A) (c : ℤ) :
    ∑ i ∈ Finset.range n, (if Q + i • P = R + P then c else 0)
      = ∑ i ∈ Finset.range n, (if Q + i • P = R then c else 0) := by
  refine (Finset.sum_nbij' (fun i => (i + 1) % n) (fun i => (i + (n - 1)) % n)
    ?_ ?_ ?_ ?_ ?_).symm
  · intro a _; exact Finset.mem_range.2 (Nat.mod_lt _ (by omega))
  · intro a _; exact Finset.mem_range.2 (Nat.mod_lt _ (by omega))
  · intro a ha; exact succ_mod_roundtrip hn (Finset.mem_range.1 ha)
  · intro a ha; exact pred_mod_roundtrip hn (Finset.mem_range.1 ha)
  · intro a _
    have hkey : Q + ((a + 1) % n) • P = (Q + a • P) + P := by
      rw [nsmul_mod n hP (a + 1), add_smul, one_smul, ← add_assoc]
    rw [hkey]
    congr 1
    simp only [eq_iff_iff]
    exact ⟨fun h => by rw [h], fun h => add_right_cancel h⟩

/-! ## ★`count` の道具 -/

section CountTools

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- `count` は有限積を和にする。 -/
theorem count_prod' {ι : Type} (v : HeightOneSpectrum W.CoordinateRing) (s : Finset ι)
    (f : ι → W.FunctionField) (hf : ∀ i ∈ s, f i ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (∏ i ∈ s, f i))
      = ∑ i ∈ s, FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (f i)) := by
  classical
  induction s using Finset.induction with
  | empty => simp [FractionalIdeal.count_one]
  | insert a s ha ih =>
      rw [Finset.prod_insert ha, Finset.sum_insert ha,
        count_mul' W v (hf a (Finset.mem_insert_self a s))
          (Finset.prod_ne_zero_iff.2 (fun i hi => hf i (Finset.mem_insert_of_mem hi))),
        ih (fun i hi => hf i (Finset.mem_insert_of_mem hi))]

/-- ★★★★★★**指数が全部 0 なら定数**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`spanSingleton z = 1` から `z` は単元の像、第 128 で定数。 -/
theorem const_of_count_eq_zero {z : W.FunctionField} (hz : z ≠ 0)
    (h : ∀ v : HeightOneSpectrum W.CoordinateRing,
      FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z) = 0) :
    ∃ c : F, c ≠ 0 ∧ z = algebraMap F W.FunctionField c := by
  have hone : FractionalIdeal.spanSingleton W.CoordinateRing⁰ z
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ 1 := by
    rw [FractionalIdeal.spanSingleton_one]
    refine eq_of_count_eq (FractionalIdeal.spanSingleton_ne_zero_iff.2 hz) one_ne_zero ?_
    intro v
    rw [h v, FractionalIdeal.count_one]
  obtain ⟨u, hu⟩ := FractionalIdeal.spanSingleton_eq_spanSingleton.1 hone
  rw [Units.smul_def, Algebra.smul_def] at hu
  have hzu : z = algebraMap W.CoordinateRing W.FunctionField
      ((u⁻¹ : W.CoordinateRingˣ) : W.CoordinateRing) := by
    have hprod : algebraMap W.CoordinateRing W.FunctionField
          ((u⁻¹ : W.CoordinateRingˣ) : W.CoordinateRing)
        * algebraMap W.CoordinateRing W.FunctionField
          ((u : W.CoordinateRingˣ) : W.CoordinateRing) = 1 := by
      rw [← map_mul]; simp
    calc z = 1 * z := (one_mul z).symm
    _ = (algebraMap W.CoordinateRing W.FunctionField
            ((u⁻¹ : W.CoordinateRingˣ) : W.CoordinateRing)
          * algebraMap W.CoordinateRing W.FunctionField
            ((u : W.CoordinateRingˣ) : W.CoordinateRing)) * z := by rw [hprod]
    _ = algebraMap W.CoordinateRing W.FunctionField
            ((u⁻¹ : W.CoordinateRingˣ) : W.CoordinateRing)
          * (algebraMap W.CoordinateRing W.FunctionField
            ((u : W.CoordinateRingˣ) : W.CoordinateRing) * z) := by rw [mul_assoc]
    _ = _ := by rw [hu, mul_one]
  obtain ⟨c, hc0, hcu⟩ := isUnit_coordinateRing (u⁻¹ : W.CoordinateRingˣ).isUnit
  refine ⟨c, hc0, ?_⟩
  rw [hzu, hcu, ← IsScalarTower.algebraMap_apply]

end CountTools

/-! ## ★平行移動の自己同型(`T = O` 込み) -/

section Translate

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- `τ` が点 `T` による平行移動であること(`T = O` なら恒等)。 -/
def IsTranslate (τ : W.FunctionField ≃ₐ[F] W.FunctionField) : W.Point → Prop
  | 0 => τ = AlgEquiv.refl
  | Point.some x₀ y₀ _ =>
      τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀

theorem isTranslate_zero {τ : W.FunctionField ≃ₐ[F] W.FunctionField} (h : τ = AlgEquiv.refl) :
    IsTranslate W τ 0 := h

theorem isTranslate_some {τ : W.FunctionField ≃ₐ[F] W.FunctionField} {x₀ y₀ : F}
    (hQ : W.Nonsingular x₀ y₀) (hx : τ (coordX W) = translateX W x₀ y₀)
    (hy : τ (coordY W) = translateY W x₀ y₀) :
    IsTranslate W τ (Point.some x₀ y₀ hQ) := ⟨hx, hy⟩

/-- ★★★**どの点についても平行移動の自己同型がある**。 -/
theorem exists_isTranslate [IsAlgClosed F] [Infinite F] [DecidableEq F] [W.IsElliptic]
    (h4 : (4 : F) ≠ 0) (T : W.Point) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField, IsTranslate W τ T := by
  match T with
  | 0 => exact ⟨AlgEquiv.refl, rfl⟩
  | Point.some x₀ y₀ hQ =>
      obtain ⟨τ, hx, hy⟩ := exists_translateAut_all W h4 hQ
      exact ⟨τ, hx, hy⟩

end Translate

/-! ## ★★★★★★★平行移動した `f_P` の指数 -/

section CountTranslate

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★**素点が等しいことと点が等しいことは同値**。 -/
theorem placeOf_eq_iff (h2 : IsUnit (2 : F)) {S T : W.Point} (hS : S ≠ 0) (hT : T ≠ 0) :
    placeOf W S hS = placeOf W T hT ↔ S = T := by
  constructor
  · intro hc
    have hcc := congrArg (pointOf W h2) hc
    rwa [pointOf_placeOf W h2 S hS, pointOf_placeOf W h2 T hT] at hcc
  · rintro rfl; rfl

/-- ★★★★★★**平行移動した `f_P` の指数(`T ≠ O`)**。 -/
theorem count_translate_fP (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hxτ : τ (coordX W) = translateX W x₀ y₀T) (hyτ : τ (coordY W) = translateY W x₀ y₀T)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hPt : n • (Point.some x y hP) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (τ (algebraMap W.CoordinateRing W.FunctionField fP)))
      = (if pointOf W h2 v + Point.some x₀ y₀T hQ = Point.some x y hP then (n : ℤ) else 0)
        - (if pointOf W h2 v + Point.some x₀ y₀T hQ = 0 then (n : ℤ) else 0) := by
  classical
  by_cases hz : pointOf W h2 v + Point.some x₀ y₀T hQ = 0
  · rw [if_neg (by rw [hz]; simp), if_pos hz]
    rw [count_translate_infty W h2 v (pointEq W h2 v) (pointAsIdeal W h2 v) τ hQ hxτ hyτ hz
      hP n hPt fP hfP]
    ring
  · rw [if_neg hz, sub_zero]
    obtain ⟨c', y₀', hns', hsum⟩ := exists_some_of_ne_zero W _ hz
    have hne' : Point.some c' y₀' hns' ≠ 0 := by simp
    set v' := placeOf W (Point.some c' y₀' hns') hne' with hv'def
    have hv' : v'.asIdeal = CoordinateRing.XYIdeal W c' (Polynomial.C y₀') :=
      placeOf_some W hns' hne'
    have hfP0 : fP ≠ 0 := fP_ne_zero W n fP hfP
    have hz0 : algebraMap W.CoordinateRing W.FunctionField fP ≠ 0 := fun h0 =>
      hfP0 (IsFractionRing.injective W.CoordinateRing W.FunctionField (by rw [h0, map_zero]))
    have hneP : Point.some x y hP ≠ 0 := by simp
    have hps : pointSpectrum W hP.1 = placeOf W (Point.some x y hP) hneP := rfl
    have hiff : (pointSpectrum W hP.1 = v')
        ↔ (Point.some x y hP = Point.some c' y₀' hns') := by
      rw [hps, hv'def]; exact placeOf_eq_iff W h2 hneP hne'
    rw [count_translate_gen W h2 v v' (pointEq W h2 v) (pointAsIdeal W h2 v) hns'.1 hv' τ hQ
      hxτ hyτ (by rw [← hsum]; rfl) hz0, count_fP W v' hP n fP hfP, hsum]
    by_cases heq : Point.some x y hP = Point.some c' y₀' hns'
    · rw [if_pos (hiff.2 heq), if_pos heq.symm]
    · rw [if_neg (fun hc => heq (hiff.1 hc)), if_neg (fun hc => heq hc.symm)]

/-- ★★★★★★★**平行移動した `f_P` の指数(全素点・`T = O` も込み)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`T = O` のときは `τ = id` で、第 2 項は `Q_v ≠ O` により 0。 -/
theorem count_translate_fP_gen (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (T : W.Point) (hT : IsTranslate W τ T)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hPt : n • (Point.some x y hP) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (τ (algebraMap W.CoordinateRing W.FunctionField fP)))
      = (if pointOf W h2 v + T = Point.some x y hP then (n : ℤ) else 0)
        - (if pointOf W h2 v + T = 0 then (n : ℤ) else 0) := by
  classical
  match T, hT with
  | 0, hT =>
      simp only [add_zero]
      rw [show τ = AlgEquiv.refl from hT, if_neg (pointOf_ne_zero W h2 v), sub_zero]
      show FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (algebraMap W.CoordinateRing W.FunctionField fP)) = _
      rw [count_fP W v hP n fP hfP]
      have hneP : Point.some x y hP ≠ 0 := by simp
      have hps : pointSpectrum W hP.1 = placeOf W (Point.some x y hP) hneP := rfl
      by_cases heq : pointOf W h2 v = Point.some x y hP
      · rw [if_pos heq, if_pos]
        rw [hps]
        exact ((placeOf_eq_iff W h2 hneP (pointOf_ne_zero W h2 v)).2 heq.symm).trans
          (placeOf_pointOf W h2 v)
      · rw [if_neg heq, if_neg]
        intro hcon
        rw [hps] at hcon
        have hcc := congrArg (pointOf W h2) hcon
        rw [pointOf_placeOf W h2 _ hneP] at hcc
        exact heq hcc.symm
  | Point.some x₀ y₀T hQ, hT =>
      exact count_translate_fP W h2 v τ hQ hT.1 hT.2 hP n hPt fP hfP

end CountTranslate

/-! ## ★★★★★★★★積は定数 -/

section Prod

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-- 点 `T` による平行移動の自己同型(選択)。 -/
noncomputable def translateAut (h4 : (4 : F) ≠ 0) (T : W.Point) :
    W.FunctionField ≃ₐ[F] W.FunctionField :=
  (exists_isTranslate W h4 T).choose

theorem isTranslate_translateAut (h4 : (4 : F) ≠ 0) (T : W.Point) :
    IsTranslate W (translateAut W h4 T) T :=
  (exists_isTranslate W h4 T).choose_spec

/-- ★★★★★★★**`∏_{i<n} τ_{iP}^*(f_P)`**。 -/
noncomputable def prodTranslate (h4 : (4 : F) ≠ 0) (P : W.Point) (n : ℕ)
    (fP : W.CoordinateRing) : W.FunctionField :=
  ∏ i ∈ Finset.range n,
    translateAut W h4 (i • P) (algebraMap W.CoordinateRing W.FunctionField fP)

theorem prodTranslate_ne_zero (h4 : (4 : F) ≠ 0) (P : W.Point) (n : ℕ)
    {fP : W.CoordinateRing} (hfP : fP ≠ 0) : prodTranslate W h4 P n fP ≠ 0 := by
  have hz : algebraMap W.CoordinateRing W.FunctionField fP ≠ 0 := fun h0 =>
    hfP (IsFractionRing.injective W.CoordinateRing W.FunctionField (by rw [h0, map_zero]))
  refine Finset.prod_ne_zero_iff.2 (fun i _ => ?_)
  intro h0
  exact hz ((translateAut W h4 (i • P)).injective (by rw [h0, map_zero]))

/-- ★★★★★★★★**積の指数はどの素点でも 0**。 -/
theorem count_prodTranslate (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (v : HeightOneSpectrum W.CoordinateRing)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hPt : n • (Point.some x y hP) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (prodTranslate W h4 (Point.some x y hP) n fP)) = 0 := by
  classical
  have hfP0 : fP ≠ 0 := fP_ne_zero W n fP hfP
  have hz : algebraMap W.CoordinateRing W.FunctionField fP ≠ 0 := fun h0 =>
    hfP0 (IsFractionRing.injective W.CoordinateRing W.FunctionField (by rw [h0, map_zero]))
  rw [prodTranslate, count_prod' W v _ _ (fun i _ h0 =>
    hz ((translateAut W h4 (i • Point.some x y hP)).injective (by rw [h0, map_zero])))]
  have hterm : ∀ i ∈ Finset.range n, FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (translateAut W h4 (i • Point.some x y hP)
          (algebraMap W.CoordinateRing W.FunctionField fP)))
      = (if pointOf W h2 v + i • Point.some x y hP = Point.some x y hP then (n : ℤ) else 0)
        - (if pointOf W h2 v + i • Point.some x y hP = 0 then (n : ℤ) else 0) := by
    intro i _
    exact count_translate_fP_gen W h2 v _ _ (isTranslate_translateAut W h4 _) hP n hPt fP hfP
  rw [Finset.sum_congr rfl hterm, Finset.sum_sub_distrib]
  have hshift := sum_shift_eq (A := W.Point) n hn hPt (pointOf W h2 v) 0 (n : ℤ)
  rw [zero_add] at hshift
  rw [hshift, sub_self]

/-- ★★★★★★★★**`∏_{i<n} τ_{iP}^*(f_P)` は定数**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★交代性 `e_n(P,P) = 1` の中核である。 -/
theorem prodTranslate_const (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hPt : n • (Point.some x y hP) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ c : F, c ≠ 0 ∧ prodTranslate W h4 (Point.some x y hP) n fP
      = algebraMap F W.FunctionField c :=
  const_of_count_eq_zero W (prodTranslate_ne_zero W h4 _ n (fP_ne_zero W n fP hfP))
    (fun v => count_prodTranslate W h2 h4 v hP n hn hPt fP hfP)

end Prod

/-! ## ★出典の紐付け(`.src`) -/

def count_translate_fP_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——平行移動した f_P の指数)",
    sectionId := "genell-thm-3-8" }

def prodTranslate_const.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——平行移動の積が定数であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
