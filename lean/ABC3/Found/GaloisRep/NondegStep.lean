import ABC3.Found.GaloisRep.TorsionAction

/-!
# Galois (G5) 第 197 ブロック —— **★★★★★★★★★非退化性を 1 つの入力に絞った**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★残る入力は 1 つだけ

非退化性の古典的な証明のうち、**`F(E)^{E[n]} = [n]^* F(E)` 以外は全部埋めた**:

    e_n(P, ·) ≡ 1
      ⟹ τ_Q(g) = g (∀Q ∈ E[n])          ← 第 195 の `algebraMap_weilPairingVal`
      ⟹ g ∈ [n]^* F(E)                   ← ★**これだけが残る**
      ⟹ v^n = f_P                         ← `μ̃` の単射性
      ⟹ (f_P) = (v)^n が分数イデアルで   ← 本ブロック
      ⟹ XYIdeal(P) = (v)                  ← `count` を `n` で割る
      ⟹ toClass P = 0 ⟹ P = O            ← mathlib の `toClass_eq_zero`

★本ブロックは `hfix` という仮定の形で「`E[n]` 不変なら `[n]^*` の像に入る」を受け取り、
**それ以外を全部証明する**。

### ★★★★★★`count` を `n` で割るところ

`(f_P) = (v)^n` と `(f_P) = XYIdeal(P)^n` から、各素点で
`n·count_v(XYIdeal(P)) = n·count_v(v)`。★`n ≠ 0` で割って第 175 の `eq_of_count_eq`。

### ★残る 1 行の中身(第 196 の測定)

第 196 で `[F(E) : F(E)^{E[n]}] = n²`(Artin)と `[n]^*F(E) ⊆ F(E)^{E[n]}` を出した。
★等号には **`[F(E) : [n]^*F(E)] = n²`**、すなわち `deg[n] = n²` が要る。
★★mathlib には分点多項式 `Φ_n`・`ΨSq_n` と**その次数**(`natDegree_Φ_le`・`coeff_Φ`・
`natDegree_ΨSq`)があり、Lüroth の副産物 `finrank_eq_max_natDegree`
(`K⟮X⟯/K⟮f⟯` の次数 = `max(deg num, deg denom)`)もある。
★★★足りないのは **`x([n]P) = Φ_n(x)/ΨSq_n(x)`**——分点多項式と群法則を結ぶ 1 本だけで、
mathlib には無い(2026-08-20 実測、`ΨSq`/`Φ` は `DivisionPolynomial/` の外に 0 件)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eq_zero_of_fP_isPow` | ★★★★★★★**`f_P` が `n` 乗なら `P = O`** |
| `eq_zero_of_weilPairingVal_trivial` | ★★★★★★★★★**`e_n(P,·) ≡ 1 ⟹ P = O`**(`hfix` 仮定) |
| `nondeg_of_fixedField` | ★★★★★★★★★**非退化性**(`hfix` 仮定) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★★★**`f_P` が `n` 乗なら `P = O`**——非退化性の最後の段。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`count` を `n` で割って `XYIdeal(P) = (v)`、`toClass` の単射性で `P = O`。 -/
theorem eq_zero_of_fP_isPow (n : ℕ) (hn : 1 ≤ n) {x y : F} (hP : W.Nonsingular x y)
    (fP : W.CoordinateRing) (hfP0 : fP ≠ 0)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {w : W.FunctionField}
    (hw : w ^ n = algebraMap W.CoordinateRing W.FunctionField fP) :
    Point.some x y hP = 0 := by
  have hw0 : w ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hw
    exact hfP0 (IsFractionRing.injective W.CoordinateRing W.FunctionField
      (by rw [← hw, map_zero]))
  have hcoe : ((CoordinateRing.XYIdeal W x (Polynomial.C y) :
      FractionalIdeal W.CoordinateRing⁰ W.FunctionField)) ^ n
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ w ^ n := by
    rw [← FractionalIdeal.coeIdeal_pow, hfP, FractionalIdeal.coeIdeal_span_singleton, ← hw,
      FractionalIdeal.spanSingleton_pow]
  have hIne : ((CoordinateRing.XYIdeal W x (Polynomial.C y) :
      FractionalIdeal W.CoordinateRing⁰ W.FunctionField)) ≠ 0 := by
    rw [← CoordinateRing.XYIdeal'_eq hP]
    exact Units.ne_zero _
  have hSne : FractionalIdeal.spanSingleton W.CoordinateRing⁰ w ≠ 0 :=
    FractionalIdeal.spanSingleton_ne_zero_iff.2 hw0
  have heq : ((CoordinateRing.XYIdeal W x (Polynomial.C y) :
      FractionalIdeal W.CoordinateRing⁰ W.FunctionField))
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ w := by
    refine eq_of_count_eq hIne hSne (fun v => ?_)
    have hc := congrArg (FractionalIdeal.count W.FunctionField v) hcoe
    rw [FractionalIdeal.count_pow, FractionalIdeal.count_pow] at hc
    have hnz : (n : ℤ) ≠ 0 := by exact_mod_cast (by omega : n ≠ 0)
    exact mul_left_cancel₀ hnz hc
  refine (Point.toClass_eq_zero _).1 ?_
  rw [Point.toClass_some hP]
  refine ClassGroup.mk_eq_one_iff.2 ?_
  rw [CoordinateRing.XYIdeal'_eq hP, heq]
  exact ⟨⟨w, by simp [FractionalIdeal.spanSingleton]⟩⟩

/-- ★★★★★★★★★**`e_n(P,·) ≡ 1` なら `P = O`**(`F(E)^{E[n]} = [n]^*F(E)` を仮定)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`hfix` が残る唯一の入力である(= `deg[n] = n²`、第 196 の測定を参照)。 -/
theorem eq_zero_of_weilPairingVal_trivial (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hfix : ∀ z : W.FunctionField,
      (∀ Q : W.Point, n • Q = 0 → translateAut W h4 Q z = z) → ∃ v, z = muExt W hinj v)
    {x y : F} (hPx : W.Nonsingular x y) (hPt : n • (Point.some x y hPx) = 0)
    (htriv : ∀ Q : W.Point, n • Q = 0 → weilPairingVal W n (Point.some x y hPx) Q = 1) :
    Point.some x y hPx = 0 := by
  obtain ⟨fP, hfP0, hfP⟩ := xyIdeal_pow_isPrincipal_integral hPx n hPt
  obtain ⟨g, hg⟩ := exists_nthRoot_mulByN W h2 hPx n hn hchar hPt μ hinj hμF hns hμP hμx hμy
    fP hfP
  have hfne : μ fP ≠ 0 := fun h0 => hfP0 (hinj (by rw [h0, map_zero]))
  have hg0 : g ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hg
    exact hfne hg.symm
  have hinv : ∀ Q : W.Point, n • Q = 0 → translateAut W h4 Q g = g := by
    intro Q hQ
    have hkey := algebraMap_weilPairingVal W h2 h4 n hn hinj hμF hns hμP hμx hμy hPx fP hfP0
      hfP hg hQ
    rw [htriv Q hQ, map_one] at hkey
    exact (div_eq_one_iff_eq hg0).1 hkey.symm
  obtain ⟨v, hv⟩ := hfix g hinv
  have hvn : v ^ n = algebraMap W.CoordinateRing W.FunctionField fP := by
    refine (muExt W hinj).injective ?_
    rw [map_pow, ← hv, hg, muExt_algebraMap]
  exact eq_zero_of_fP_isPow W n hn hPx fP hfP0 hfP hvn

/-- ★★★★★★★★★**非退化性**(`F(E)^{E[n]} = [n]^*F(E)` を仮定した形)。 -/
theorem nondeg_of_fixedField (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hfix : ∀ z : W.FunctionField,
      (∀ Q : W.Point, n • Q = 0 → translateAut W h4 Q z = z) → ∃ v, z = muExt W hinj v)
    {x y : F} (hPx : W.Nonsingular x y) (hPt : n • (Point.some x y hPx) = 0)
    (hne : Point.some x y hPx ≠ 0) :
    ∃ Q : W.Point, n • Q = 0 ∧ weilPairingVal W n (Point.some x y hPx) Q ≠ 1 := by
  by_contra hcon
  refine hne (eq_zero_of_weilPairingVal_trivial W h2 h4 n hn hchar hinj hμF hns hμP hμx hμy
    hfix hPx hPt (fun Q hQ => ?_))
  by_contra hval
  exact hcon ⟨Q, hQ, hval⟩

/-! ## ★出典の紐付け(`.src`) -/

def eq_zero_of_fP_isPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——f_P が n 乗なら P = O)",
    sectionId := "genell-thm-3-8" }

def nondeg_of_fixedField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
