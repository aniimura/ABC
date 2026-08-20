import ABC3.Found.GaloisRep.MuComm

/-!
# Galois (G5) 第 190 ブロック —— **★★★★★★★★★交代性 `e_n(P,P) = 1`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★5 性質の 3 つ目

    e_n(P, P) = 1

★第 186-189 で積み上げた 4 つが、ここで合流する。

| 材料 | 出どころ |
|---|---|
| `n·P' = P` なる `P'` | 第 186(位数からの割り算) |
| `∏_{i<n} τ_{iP}^*(f_P)` は定数 | 第 188(極の位数 `−n` は第 187) |
| `τ_T ∘ [n]^* = [n]^* ∘ τ_{nT}` | 第 189 |
| `定数の n 乗根は定数` | 第 177 |

### ★★★★★★機構

1. `H := ∏_{i<n} τ_{iP'}^*(g)` を作る(`g^n = [n]^* f_P`)
2. `H^n = ∏_i τ_{iP'}^*([n]^* f_P) = ∏_i [n]^*(τ_{iP}^* f_P) = [n]^*(定数) = 定数`
   ★ここで `n·(i·P') = i·(n·P') = i·P` を使う(第 189 の交換則)
3. 第 177 で `H` は定数、よって `τ_{P'}^* H = H`
4. **伸縮**: `τ_{P'}^* H = ∏_{i<n} τ_{(i+1)P'}^*(g)` なので

       (∏_{i<n} τ_{(i+1)P'}^* g) · g = H · τ_{nP'}^*(g)
       ⟹  H · g = H · τ_P^*(g)  ⟹  τ_P^*(g) = g

   ★`τ_{P'} ∘ τ_{iP'} = τ_{(i+1)P'}` は第 189 の `isTranslate_trans` と一意性から。
5. したがって `e_n(P,P) = τ_P^*(g)/g = 1`

### ★逸脱(記録)

`hchar` の範囲が `k ≤ n²` である(第 186 が `#E[n²] = n⁴` を使うため)。
★最終消費者 `det_cyclotomic` は `[CharZero K]` の下なので無害。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isTranslate_unique` | ★★★★平行移動の自己同型は一意 |
| `translateAut_trans` | ★★★合成は「点の和」の平行移動 |
| `prodRoot_pow` | ★★★★★★★★`H^n` は定数 |
| `translate_root_fixed` | ★★★★★★★★★**`τ_P^*(g) = g`** |
| `weilPairingVal_self` | ★★★★★★★★★**交代性 `e_n(P,P) = 1`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★平行移動の一意性と合成 -/

section Unique

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- ★★★★**平行移動の自己同型は一意**。 -/
theorem isTranslate_unique {τ τ' : W.FunctionField ≃ₐ[F] W.FunctionField} {T : W.Point}
    (hT : IsTranslate W τ T) (hT' : IsTranslate W τ' T) : τ = τ' := by
  match T, hT, hT' with
  | 0, hT, hT' => rw [show τ = AlgEquiv.refl from hT, show τ' = AlgEquiv.refl from hT']
  | Point.some x₀ y₀ hQ, hT, hT' =>
      exact aut_ext W (hT.1.trans hT'.1.symm) (hT.2.trans hT'.2.symm)

end Unique

section Trans

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

theorem translateAut_zero (h4 : (4 : F) ≠ 0) :
    translateAut W h4 (0 : W.Point) = AlgEquiv.refl :=
  isTranslate_translateAut W h4 0

/-- ★★★平行移動の合成は「点の和」の平行移動。 -/
theorem translateAut_trans (h4 : (4 : F) ≠ 0) (T T' : W.Point) :
    (translateAut W h4 T').trans (translateAut W h4 T) = translateAut W h4 (T + T') :=
  isTranslate_unique W (isTranslate_trans W (isTranslate_translateAut W h4 T)
    (isTranslate_translateAut W h4 T')) (isTranslate_translateAut W h4 (T + T'))

/-! ## ★★★★★★★★`H^n` は定数 -/

/-- ★★★★★★★★`H := ∏_{i<n} τ_{iP'}(g)` の `n` 乗は定数。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 189 の交換則で `[n]^*` の外に出し、第 188 の「積は定数」を当てる。 -/
theorem prodRoot_pow (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hPt : n • (Point.some x y hP) = 0)
    {P' : W.Point} (hP' : n • P' = Point.some x y hP)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {g : W.FunctionField} (hg : g ^ n = μ fP) :
    ∃ c : F, c ≠ 0 ∧ (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g) ^ n
      = algebraMap F W.FunctionField c := by
  obtain ⟨c, hc0, hconst⟩ := prodTranslate_const W h2 h4 hP n hn hPt fP hfP
  refine ⟨c, hc0, ?_⟩
  rw [← Finset.prod_pow]
  have hterm : ∀ i ∈ Finset.range n,
      (translateAut W h4 (i • P') g) ^ n
        = muExt W hinj (translateAut W h4 (i • Point.some x y hP)
            (algebraMap W.CoordinateRing W.FunctionField fP)) := by
    intro i _
    rw [← map_pow, hg, ← muExt_algebraMap W hinj fP]
    refine aut_comp_muExt_gen W hinj hμF hns hμx hμy n hμP _ _ (i • P')
      (i • Point.some x y hP) (isTranslate_translateAut W h4 _)
      (isTranslate_translateAut W h4 _) ?_ _
    rw [← hP', ← mul_smul, ← mul_smul, Nat.mul_comm]
  rw [Finset.prod_congr rfl hterm, ← map_prod, ← prodTranslate, hconst,
    muExt_const W hinj hμF c]

/-! ## ★★★★★★★★★伸縮 -/

/-- ★★★★★★★★★**`τ_P^*(g) = g`**——交代性の中身。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`H` が定数なので `τ_{P'}^* H = H`。★★左辺を伸縮すると `H·g = H·τ_P^*(g)`。 -/
theorem translate_root_fixed (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hPt : n • (Point.some x y hP) = 0)
    {P' : W.Point} (hP' : n • P' = Point.some x y hP)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {g : W.FunctionField} (hg : g ^ n = μ fP) :
    translateAut W h4 (Point.some x y hP) g = g := by
  obtain ⟨c, hc0, hpow⟩ := prodRoot_pow W h2 h4 hP n hn hPt hP' hinj hμF hns hμP hμx hμy fP hfP hg
  obtain ⟨ζ, hζ0, hζ⟩ := const_of_pow_eq_const W h2 hn hc0 hpow
  have hH0 : (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g) ≠ 0 := by
    rw [hζ]
    exact fun h0 => hζ0 ((algebraMap F W.FunctionField).injective (by rw [h0, map_zero]))
  have hfix : translateAut W h4 P' (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g)
      = ∏ i ∈ Finset.range n, translateAut W h4 (i • P') g := by
    rw [hζ]; exact (translateAut W h4 P').commutes ζ
  have hshift : ∀ i : ℕ, translateAut W h4 P' (translateAut W h4 (i • P') g)
      = translateAut W h4 ((i + 1) • P') g := by
    intro i
    have hcomp : (translateAut W h4 (i • P')).trans (translateAut W h4 P')
        = translateAut W h4 ((i + 1) • P') := by
      rw [translateAut_trans]
      congr 1
      rw [succ_nsmul, add_comm]
    exact congrFun (congrArg (fun e : W.FunctionField ≃ₐ[F] W.FunctionField => e.toFun) hcomp) g
  have hprodshift : (∏ i ∈ Finset.range n, translateAut W h4 ((i + 1) • P') g)
      = ∏ i ∈ Finset.range n, translateAut W h4 (i • P') g := by
    calc (∏ i ∈ Finset.range n, translateAut W h4 ((i + 1) • P') g)
        = ∏ i ∈ Finset.range n, translateAut W h4 P' (translateAut W h4 (i • P') g) :=
          Finset.prod_congr rfl (fun i _ => (hshift i).symm)
      _ = translateAut W h4 P' (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g) :=
          (map_prod _ _ _).symm
      _ = _ := hfix
  have h1 : ∏ k ∈ Finset.range (n + 1), translateAut W h4 (k • P') g
      = (∏ k ∈ Finset.range n, translateAut W h4 ((k + 1) • P') g)
        * translateAut W h4 ((0 : ℕ) • P') g :=
    Finset.prod_range_succ' (fun i => translateAut W h4 (i • P') g) n
  have h2' : ∏ k ∈ Finset.range (n + 1), translateAut W h4 (k • P') g
      = (∏ k ∈ Finset.range n, translateAut W h4 (k • P') g)
        * translateAut W h4 (n • P') g :=
    Finset.prod_range_succ (fun i => translateAut W h4 (i • P') g) n
  rw [hprodshift] at h1
  have hkey : (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g)
      * translateAut W h4 ((0 : ℕ) • P') g
      = (∏ i ∈ Finset.range n, translateAut W h4 (i • P') g)
      * translateAut W h4 (n • P') g := by rw [← h1, ← h2']
  rw [zero_smul, translateAut_zero W h4, hP'] at hkey
  exact (mul_left_cancel₀ hH0 hkey).symm

/-! ## ★★★★★★★★★交代性 -/

/-- ★★★★★★★★★**交代性 `e_n(P,P) = 1`**——Weil 対の 5 性質の 3 つ目。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`hchar` の範囲は `k ≤ n²`(第 186 が `#E[n²] = n⁴` を使うため)。
★★最終消費者 `det_cyclotomic` は `[CharZero K]` の下なので無害。 -/
theorem weilPairingVal_self (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n ^ 2 → (k : F) ≠ 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {x y : F} (hP : W.Nonsingular x y) (hPt : n • (Point.some x y hP) = 0) :
    weilPairingVal W n (Point.some x y hP) (Point.some x y hP) = 1 := by
  have hcharn : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0 :=
    fun k h1 hk => hchar k h1 (le_trans hk (Nat.le_self_pow two_ne_zero n))
  have hinj := mulByN_injective W h2 n hn (hcharn n hn le_rfl) μ hμF hns hμx hμy hμP
  obtain ⟨fP, hfP0, hfP⟩ := xyIdeal_pow_isPrincipal_integral hP n hPt
  obtain ⟨g, hg⟩ := exists_nthRoot_mulByN W h2 hP n hn hcharn hPt μ hinj hμF hns hμP hμx hμy
    fP hfP
  obtain ⟨P', hP'⟩ := exists_nsmul_eq_point W n hn hchar hPt
  have hfixg := translate_root_fixed W h2 h4 hP n hn hPt hP' hinj hμF hns hμP hμx hμy fP hfP hg
  have hg0 : g ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hg
    exact hfP0 (hinj (by rw [← hg, map_zero]))
  have hτ := isTranslate_translateAut W h4 (Point.some x y hP)
  refine weilPairingVal_eq W h2 hn ⟨x, y, hP, x, y, hP, fP, μ, g,
    translateAut W h4 (Point.some x y hP), xn, yn, hns, rfl, rfl, hfP, hinj, hμF,
    hμx, hμy, hμP, hg, hτ.1, hτ.2, ?_⟩
  rw [hfixg, div_self hg0, map_one]

end Trans

/-! ## ★出典の紐付け(`.src`) -/

def translate_root_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——τ_P^*(g) = g)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_self.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
