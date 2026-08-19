import ABC3.Found.GaloisRep.GoodBase

/-!
# Galois (G1) 第 65 ブロック —— **★★★★★★★`#E[n] = n²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★数え上げが閉じる

良い基点 `c`(第 64)と `Q = (c, d)` に対し

    [n]⁻¹(Q)  ⟷  g_c の相異なる根        (x 座標を取る写像)

が**全単射**である:

| 向き | 根拠 |
|---|---|
| 単射 | ★`Ψ₂Sq(c) ≠ 0` ⟹ `Q ≠ −Q`——`R` と `−R` の片方しかファイバーに入らない |
| 定義可能 | ★★条件 (3)——ファイバーに**位数 ≤ n の点が入らない** |
| 全射 | ★★★条件 (2) + 乗法公式——各根から点が作れる |

★したがって `#[n]⁻¹(Q) = n²`(第 63)。
★★平行移動 `P ↦ P + R₀` で `E[n] ≅ [n]⁻¹(Q)` なので **`#E[n] = n²`** ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `order_eq_min_root` | ★★最小の根の添字は位数 |
| `fiber_nondegenerate` | ★★★★ファイバーの点は退化しない |
| `fiber_isRoot` / `root_fiber` | ★★★両向き |
| `fiber_inj` / `fiber_image` | ★★★単射・像 |
| `fiber_ncard` | ★★★★★★`#[n]⁻¹(Q) = n²` |
| `torsion_ncard` | ★★★★★★★**`#E[n] = n²`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★点の x 座標(無限遠点は `0`)。 -/
def xCoord : W.toAffine.Point → F
  | .zero => 0
  | .some x _ _ => x

/-- ★★**最小の根の添字は点の位数である**。 -/
theorem order_eq_min_root {x y : F} (h : W.toAffine.Nonsingular x y) (K : ℕ) (hK : 1 ≤ K)
    (hroot : (W.ΨSq (K : ℤ)).eval x = 0)
    (hmin : ∀ k : ℕ, 1 ≤ k → k < K → (W.ΨSq (k : ℤ)).eval x ≠ 0) :
    addOrderOf (Point.some x y h) = K := by
  have hKz : K • (Point.some x y h) = 0 := smul_eq_zero_of_min_root' W h K hK hroot hmin
  have hdK : addOrderOf (Point.some x y h) ∣ K := addOrderOf_dvd_of_nsmul_eq_zero hKz
  have hd1 : 1 ≤ addOrderOf (Point.some x y h) :=
    Nat.pos_of_ne_zero (fun h0 => by rw [h0, Nat.zero_dvd] at hdK; omega)
  have hdle : addOrderOf (Point.some x y h) ≤ K := Nat.le_of_dvd (by omega) hdK
  by_contra hne
  have hlt : addOrderOf (Point.some x y h) < K := lt_of_le_of_ne hdle hne
  obtain ⟨x1, y1, h1', hP1, -, -⟩ := mulOK_of_ne W h (addOrderOf (Point.some x y h)) hd1
    (fun k hk1 hk2 => hmin k hk1 (by omega))
  rw [addOrderOf_nsmul_eq_zero] at hP1
  exact Point.some_ne_zero h1' hP1.symm

/-- ★★★★**ファイバーの点は退化しない**。 -/
theorem fiber_nondegenerate (n : ℕ) (hn : 2 ≤ n) (c : F)
    (hc3 : ∀ z : F, (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval z = 0) →
      ∀ j : ℕ, j < n → c ≠ (W.Φ (j : ℤ)).eval z / (W.ΨSq (j : ℤ)).eval z)
    {x y : F} (h : W.toAffine.Nonsingular x y)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ)
    (hR : n • (Point.some x y h) = Point.some c dQ hQ) :
    ∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval x ≠ 0 := by
  classical
  by_contra hcon
  push_neg at hcon
  obtain ⟨k0, hk01, hk0n, hk0z⟩ := hcon
  have hex : ∃ k : ℕ, 1 ≤ k ∧ (W.ΨSq (k : ℤ)).eval x = 0 := ⟨k0, hk01, hk0z⟩
  obtain ⟨hK1, hKroot⟩ := Nat.find_spec hex
  set K := Nat.find hex with hKdef
  have hKn : K ≤ n := le_trans (Nat.find_le ⟨hk01, hk0z⟩) hk0n
  have hmin : ∀ k : ℕ, 1 ≤ k → k < K → (W.ΨSq (k : ℤ)).eval x ≠ 0 :=
    fun k hk1 hkK hkz => Nat.find_min hex hkK ⟨hk1, hkz⟩
  have hord : addOrderOf (Point.some x y h) = K := order_eq_min_root W h K hK1 hKroot hmin
  have hKz : K • (Point.some x y h) = 0 := by rw [← hord]; exact addOrderOf_nsmul_eq_zero _
  set j := n % K with hj
  have hjR : n • (Point.some x y h) = j • (Point.some x y h) := by
    have hsplit : n = (n / K) * K + j := by rw [hj]; exact (Nat.div_add_mod' n K).symm
    conv_lhs => rw [hsplit]
    rw [add_smul, mul_smul, hKz, smul_zero, zero_add]
  have hj0 : j ≠ 0 := by
    intro h0
    rw [h0, zero_smul] at hjR
    rw [hjR] at hR
    exact Point.some_ne_zero hQ hR.symm
  have hjK : j < K := Nat.mod_lt _ (by omega)
  have hjne : (W.ΨSq (j : ℤ)).eval x ≠ 0 := hmin j (by omega) hjK
  obtain ⟨x1, y1, h1', hP1, hx1, -⟩ := mulOK_of_ne W h j (by omega)
    (fun k hk1 hk2 => hmin k hk1 (by omega))
  rw [← hjR, hR] at hP1
  obtain ⟨hxx, -⟩ := Point.some.inj hP1.symm
  refine hc3 x ⟨K, hK1, hKn, hKroot⟩ j (by omega) ?_
  rw [← hxx]
  field_simp
  linear_combination hx1

/-- ★★★ファイバーの点の x 座標は `g_c` の根。 -/
theorem fiber_isRoot (n : ℕ) (hn : 2 ≤ n) (c : F)
    (hc3 : ∀ z : F, (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval z = 0) →
      ∀ j : ℕ, j < n → c ≠ (W.Φ (j : ℤ)).eval z / (W.ΨSq (j : ℤ)).eval z)
    {x y : F} (h : W.toAffine.Nonsingular x y)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ)
    (hR : n • (Point.some x y h) = Point.some c dQ hQ) :
    (gPoly W (n : ℤ) c).IsRoot x := by
  have hnd := fiber_nondegenerate W n hn c hc3 h hQ hR
  obtain ⟨x', y', h'', hP, hx', -⟩ := mulOK_of_ne W h n (by omega) hnd
  rw [hR] at hP
  obtain ⟨rfl, -⟩ := Point.some.inj hP.symm
  simp only [gPoly, Polynomial.IsRoot, eval_sub, eval_mul, eval_C, sub_eq_zero]
  linear_combination -hx'

/-- ★★★`g_c` の根は必ずファイバーの点から来る。 -/
theorem root_fiber (n : ℕ) (hn : 2 ≤ n) (c r : F)
    (hgood : ∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval r ≠ 0)
    (hroot : (gPoly W (n : ℤ) c).IsRoot r)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ)
    {y0 : F} (h0 : W.toAffine.Nonsingular r y0) :
    ∃ (y : F) (h : W.toAffine.Nonsingular r y),
      n • (Point.some r y h) = Point.some c dQ hQ := by
  obtain ⟨x', y', h'', hP, hx', -⟩ := mulOK_of_ne W h0 n (by omega) hgood
  have hcc : x' = c := by
    have hh : (W.Φ (n : ℤ)).eval r = c * (W.ΨSq (n : ℤ)).eval r := by
      have hrr := hroot
      simp only [gPoly, Polynomial.IsRoot, eval_sub, eval_mul, eval_C, sub_eq_zero] at hrr
      exact hrr
    refine mul_right_cancel₀ (hgood n (by omega) le_rfl) ?_
    rw [hx', hh]
  rcases (Point.X_eq_iff (h₁ := h'') (h₂ := hQ)).mp hcc with hcase | hcase
  · exact ⟨y0, h0, by rw [hP, hcase]⟩
  · refine ⟨W.toAffine.negY r y0, (nonsingular_neg ..).mpr h0, ?_⟩
    rw [← Point.neg_some h0, smul_neg, hP, hcase, neg_neg]

/-! ## ★★★★★★数え上げ -/

theorem fiber_inj (n : ℕ) (c : F) (hc1 : W.Ψ₂Sq.eval c ≠ 0)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ) :
    Set.InjOn (xCoord W) {R : W.toAffine.Point | n • R = Point.some c dQ hQ} := by
  have hQne : dQ ≠ W.toAffine.negY c dQ := fun hh =>
    hc1 ((psi2Sq_eval_eq_zero_iff W hQ.left).2 hh)
  rintro (_ | ⟨x1, y1, h1⟩) hR1 (_ | ⟨x2, y2, h2⟩) hR2 hxx
  · rfl
  · exfalso
    have h0 : n • (0 : W.toAffine.Point) = Point.some c dQ hQ := hR1
    rw [smul_zero] at h0
    exact Point.some_ne_zero hQ h0.symm
  · exfalso
    have h0 : n • (0 : W.toAffine.Point) = Point.some c dQ hQ := hR2
    rw [smul_zero] at h0
    exact Point.some_ne_zero hQ h0.symm
  · simp only [xCoord] at hxx
    rcases (Point.X_eq_iff (h₁ := h1) (h₂ := h2)).mp hxx with hcase | hcase
    · exact hcase
    · exfalso
      have e1 : n • (Point.some x1 y1 h1) = Point.some c dQ hQ := hR1
      have e2 : n • (Point.some x2 y2 h2) = Point.some c dQ hQ := hR2
      rw [hcase, smul_neg, e2, Point.neg_some] at e1
      obtain ⟨-, hy⟩ := Point.some.inj e1
      exact hQne hy.symm

theorem fiber_image [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 2 ≤ n) (c : F)
    (hc2 : ∀ r : F, (gPoly W (n : ℤ) c).IsRoot r →
      (¬ (derivative (gPoly W (n : ℤ) c)).IsRoot r)
      ∧ (∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval r ≠ 0))
    (hc3 : ∀ z : F, (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval z = 0) →
      ∀ j : ℕ, j < n → c ≠ (W.Φ (j : ℤ)).eval z / (W.ΨSq (j : ℤ)).eval z)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ) :
    (xCoord W) '' {R : W.toAffine.Point | n • R = Point.some c dQ hQ}
      = ↑((gPoly W (n : ℤ) c).roots.toFinset) := by
  classical
  have hnabs : 2 ≤ (n : ℤ).natAbs := by simpa using hn
  have hgne : gPoly W (n : ℤ) c ≠ 0 := (gPoly_monic W (n : ℤ) hnabs c).ne_zero
  ext r
  constructor
  · rintro ⟨R, hR, rfl⟩
    rcases R with _ | ⟨x1, y1, h1⟩
    · exfalso
      have h0 : n • (0 : W.toAffine.Point) = Point.some c dQ hQ := hR
      rw [smul_zero] at h0
      exact Point.some_ne_zero hQ h0.symm
    · simp only [xCoord, Finset.mem_coe, Multiset.mem_toFinset, Polynomial.mem_roots hgne]
      exact fiber_isRoot W n hn c hc3 h1 hQ hR
  · intro hr
    simp only [Finset.mem_coe, Multiset.mem_toFinset, Polynomial.mem_roots hgne] at hr
    obtain ⟨y0, h0⟩ := exists_nonsingular W hΔ r
    obtain ⟨y, hy, hfib⟩ := root_fiber W n hn c r (hc2 r hr).2 hr hQ h0
    exact ⟨Point.some r y hy, hfib, rfl⟩

/-- ★★★★★★**ファイバーの個数は `n²`**。 -/
theorem fiber_ncard [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 2 ≤ n) (c : F)
    (hc1 : W.Ψ₂Sq.eval c ≠ 0)
    (hc2 : ∀ r : F, (gPoly W (n : ℤ) c).IsRoot r →
      (¬ (derivative (gPoly W (n : ℤ) c)).IsRoot r)
      ∧ (∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval r ≠ 0))
    (hc3 : ∀ z : F, (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval z = 0) →
      ∀ j : ℕ, j < n → c ≠ (W.Φ (j : ℤ)).eval z / (W.ΨSq (j : ℤ)).eval z)
    {dQ : F} (hQ : W.toAffine.Nonsingular c dQ) :
    {R : W.toAffine.Point | n • R = Point.some c dQ hQ}.ncard = n ^ 2 := by
  classical
  have hnabs : 2 ≤ (n : ℤ).natAbs := by simpa using hn
  have hcard := gPoly_toFinset_card W (n : ℤ) hnabs c (fun r hr => (hc2 r hr).1)
  rw [show (n : ℤ).natAbs = n by simp] at hcard
  rw [← hcard, ← Set.ncard_coe_finset, ← fiber_image W hΔ n hn c hc2 hc3 hQ,
    Set.InjOn.ncard_image (fiber_inj W n c hc1 hQ)]

/-- ★★★★★★★**`#E[n] = n²`**(代数閉体、`k ≤ n` で `(k : F) ≠ 0`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが `E[n] ≅ (ℤ/n)²` の数の側である。 -/
theorem torsion_ncard [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 2 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    {P : W.toAffine.Point | n • P = 0}.ncard = n ^ 2 := by
  classical
  obtain ⟨c, hc1, hc2, hc3⟩ := exists_good_c' W hΔ n hn hchar
  obtain ⟨dQ, hQ⟩ := exists_nonsingular W hΔ c
  have hfib := fiber_ncard W hΔ n hn c hc1 hc2 hc3 hQ
  have hne : {R : W.toAffine.Point | n • R = Point.some c dQ hQ}.Nonempty :=
    Set.nonempty_of_ncard_ne_zero (by rw [hfib]; positivity)
  obtain ⟨R0, hR0⟩ := hne
  have himg : (fun P => P + R0) '' {P : W.toAffine.Point | n • P = 0}
      = {R : W.toAffine.Point | n • R = Point.some c dQ hQ} := by
    ext R
    constructor
    · rintro ⟨P, hP, rfl⟩
      show n • (P + R0) = _
      rw [smul_add, (by exact hP : n • P = 0), zero_add]
      exact hR0
    · intro hR
      refine ⟨R - R0, ?_, sub_add_cancel R R0⟩
      show n • (R - R0) = 0
      rw [smul_sub, (by exact hR : n • R = _), (by exact hR0 : n • R0 = _), sub_self]
  rw [← hfib, ← himg, Set.ncard_image_of_injective _ (fun a b hab => by
    simpa using hab : Function.Injective (fun P : W.toAffine.Point => P + R0))]

/-! ## ★出典の紐付け(`.src`) -/

def torsion_ncard.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れ部分群の位数が n²)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
