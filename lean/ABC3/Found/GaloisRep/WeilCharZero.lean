import ABC3.Found.GaloisRep.WeilAlt

/-!
# Galois (G5) 第 191 ブロック —— **★★★★★★★★3 性質を標数 0 で完全形にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★スケルトンに渡せる形にする

第 179・184・190 の 3 性質は、点が `Point.some` であることや `μ` のデータを
**仮定に持っている**。★スケルトン `Skeleton/GaloisRep/WeilPairingDef.lean` が要求するのは
`P Q : W.Point` の**任意の**点についての形である。

### ★★★標数 0 なら余計な仮定は全部消える

最終消費者 `det_galRep_eq_cyclotomic` は `[CharZero K]`・`[IsAlgClosed L]` の下にある。
★`[CharZero F]` から `IsUnit (2:F)`・`(4:F) ≠ 0`・`hchar`・`Infinite F` がすべて出る。
★★`μ` は第 125 の `exists_mulByNHom_alg_charZero` で**仮定なしに存在する**。

### ★★★★★★★`O` の場合と `P₂ = −P₁` の場合

| 場合 | 扱い |
|---|---|
| `P = O` または `Q = O` | `WeilSpec` は `P = Point.some …` を要求するので witness が無い ⟹ 値は `1` |
| `P₁ + P₂ = O` | 第 184 が使えない(和が `Point.some` でない)——**別に証明する** |

★`P₁ + P₂ = O` の場合は `f_P · f_{−P} = c·(x − x_P)^n`(第 141)から出る:

    (g₁g₂/z)^n = c   (z := μ(x) − x_P)  ⟹  g₁g₂ = ζ·z

★★`τ_Q` は `μ` の像を固定する(`n·Q = O`、第 168)ので **`τ_Q z = z`**、
よって第 181 の `aut_div_mul` で `(τg₁/g₁)(τg₂/g₂) = 1`。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `weilPairingVal_zero_left` / `_right` | `O` を入れると値は 1 |
| `weilPairingVal_mul_neg` | ★★★★★★★**`e_n(P,Q)·e_n(−P,Q) = 1`** |
| `weilPairingVal_add_left` | ★★★★★★★★**双線型性(完全形)** |
| `weilPairingVal_alt` | ★★★★★★★★★**交代性(完全形)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★`O` を入れたときの値 -/

section Zero

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

theorem weilPairingVal_zero_left (n : ℕ) (Q : W.Point) : weilPairingVal W n 0 Q = 1 := by
  refine weilPairingVal_of_not W n 0 Q ?_
  rintro ⟨c, x, y, hP, x₀, y₀, hQ, fP, μ, g, τ, xn, yn, hnsg, hPeq, -⟩
  exact absurd hPeq (by simp)

theorem weilPairingVal_zero_right (n : ℕ) (P : W.Point) : weilPairingVal W n P 0 = 1 := by
  refine weilPairingVal_of_not W n P 0 ?_
  rintro ⟨c, x, y, hP, x₀, y₀, hQ, fP, μ, g, τ, xn, yn, hnsg, hPeq, hQeq, -⟩
  exact absurd hQeq (by simp)

end Zero

/-! ## ★★★★★★★`e_n(P,Q)·e_n(−P,Q) = 1` -/

section Neg

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★★★**`e_n(P,Q)·e_n(−P,Q) = 1`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 141 の `f_P · f_{−P} = c(x − x_P)^n` と、`τ_Q` が `μ` の像を固定すること(第 168)。 -/
theorem weilPairingVal_mul_neg (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {x y : F} (hP : W.Nonsingular x y) (hPt : n • (Point.some x y hP) = 0)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (hQt : n • (Point.some x₀ y₀ hQ) = 0) :
    weilPairingVal W n (Point.some x y hP) (Point.some x₀ y₀ hQ)
      * weilPairingVal W n (Point.some x (W.negY x y) ((nonsingular_neg x y).mpr hP))
        (Point.some x₀ y₀ hQ) = 1 := by
  have hinj := mulByN_injective W h2 n hn (hchar n hn le_rfl) μ hμF hns hμx hμy hμP
  set hN := (nonsingular_neg x y).mpr hP with hNdef
  have hNt : n • (Point.some x (W.negY x y) hN) = 0 := by
    have hneg : Point.some x (W.negY x y) hN = -(Point.some x y hP) := rfl
    rw [hneg, smul_neg, hPt, neg_zero]
  obtain ⟨fP, hfP0, hfP⟩ := xyIdeal_pow_isPrincipal_integral hP n hPt
  obtain ⟨fN, hfN0, hfN⟩ := xyIdeal_pow_isPrincipal_integral hN n hNt
  obtain ⟨g₁, hg₁⟩ := exists_nthRoot_mulByN W h2 hP n hn hchar hPt μ hinj hμF hns hμP hμx hμy
    fP hfP
  obtain ⟨g₂, hg₂⟩ := exists_nthRoot_mulByN W h2 hN n hn hchar hNt μ hinj hμF hns hμP hμx hμy
    fN hfN
  obtain ⟨τ, hτx, hτy⟩ := exists_translateAut_all W h4 hQ
  have hcomp := aut_comp_mulByN W τ hQ hτx hτy n hQt μ hμF hns hμx hμy hμP
  have hτz : ∀ r : W.CoordinateRing, τ (μ r) = μ r := fun r =>
    congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hcomp) r
  have hne : ∀ f : W.CoordinateRing, f ≠ 0 → μ f ≠ 0 :=
    fun f hf h0 => hf (hinj (by rw [h0, map_zero]))
  set z := μ (genX W) - algebraMap F W.FunctionField x with hzdef
  have hz : z = μ (genX W - algebraMap F W.CoordinateRing x) := by
    rw [hzdef, map_sub, hμF x]
  have hzne : z ≠ 0 := by
    rw [hz]
    exact hne _ (fun h0 => genX_ne_const W x (by rwa [sub_eq_zero] at h0))
  have hτzz : τ z = z := by rw [hz]; exact hτz _
  obtain ⟨c, hc0, hprod⟩ := mu_fP_mul_mu_fNegP W hP n fP fN hfP hfN μ hμF
  have hpow : (g₁ * g₂ / z) ^ n = algebraMap F W.FunctionField c := by
    rw [div_pow, mul_pow, hg₁, hg₂, hprod, ← hzdef]
    field_simp
  obtain ⟨ζ, hζ0, hζ⟩ := const_of_pow_eq_const W h2 hn hc0 hpow
  have heq : g₁ * g₂ = algebraMap F W.FunctionField ζ * z * 1 := by
    have hd := (div_eq_iff hzne).1 hζ
    rw [hd, mul_one]
  have hcore := aut_div_mul W τ hζ0 hzne one_ne_zero hτzz heq
  rw [map_one, div_one] at hcore
  obtain ⟨c₁, hc₁0, hd₁⟩ := exists_const_aut_div W h2 hn hg₁ (hne fP hfP0) (hτz fP)
  obtain ⟨c₂, hc₂0, hd₂⟩ := exists_const_aut_div W h2 hn hg₂ (hne fN hfN0) (hτz fN)
  have hval : ∀ {a b : F} (hPx : W.Nonsingular a b) (f : W.CoordinateRing)
      (hf : (CoordinateRing.XYIdeal W a (Polynomial.C b)) ^ n = Ideal.span {f})
      (g : W.FunctionField) (hg : g ^ n = μ f) {d : F}
      (hd : τ g / g = algebraMap F W.FunctionField d),
      weilPairingVal W n (Point.some a b hPx) (Point.some x₀ y₀ hQ) = d := by
    intro a b hPx f hf g hg d hd
    exact weilPairingVal_eq W h2 hn
      ⟨a, b, hPx, x₀, y₀, hQ, f, μ, g, τ, xn, yn, hns, rfl, rfl, hf, hinj, hμF,
        hμx, hμy, hμP, hg, hτx, hτy, hd⟩
  rw [hval hP fP hfP g₁ hg₁ hd₁, hval hN fN hfN g₂ hg₂ hd₂]
  rw [hd₁, hd₂, ← map_mul] at hcore
  exact (algebraMap F W.FunctionField).injective (hcore.trans (map_one _).symm)

end Neg

/-! ## ★★★標数 0 で余計な仮定が消える -/

section CharZero

variable {F : Type} [Field F] [DecidableEq F] [CharZero F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

theorem isUnit_two_of_charZero : IsUnit (2 : F) := Ne.isUnit (two_ne_zero (α := F))

theorem four_ne_zero_of_charZero : (4 : F) ≠ 0 := by
  have h : ((4 : ℕ) : F) ≠ 0 := Nat.cast_ne_zero.2 (by norm_num)
  simp only [Nat.cast_ofNat] at h
  exact h

theorem natCast_ne_zero_of_charZero (n : ℕ) : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0 :=
  fun k hk _ => Nat.cast_ne_zero.2 (by omega)

/-- ★★★★★★★★**双線型性(第 1 変数、標数 0・完全版)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`O` の場合と `P₁ + P₂ = O` の場合を別扱いにして、任意の点で成り立つ形にした。 -/
theorem weilPairingVal_add_left (n : ℕ) (hn : 1 ≤ n) (P₁ P₂ Q : W.Point)
    (h₁ : n • P₁ = 0) (h₂ : n • P₂ = 0) (hQt : n • Q = 0) :
    weilPairingVal W n (P₁ + P₂) Q
      = weilPairingVal W n P₁ Q * weilPairingVal W n P₂ Q := by
  haveI := isDedekindDomain_coordinateRing (isUnit_two_of_charZero (F := F)) W
  match Q, hQt with
  | 0, _ => rw [weilPairingVal_zero_right, weilPairingVal_zero_right,
      weilPairingVal_zero_right, one_mul]
  | Point.some x₀ y₀ hQ, hQt =>
    match P₁, h₁ with
    | 0, _ => rw [zero_add, weilPairingVal_zero_left, one_mul]
    | Point.some x₁ y₁ hP₁, h₁ =>
      match P₂, h₂ with
      | 0, _ => rw [add_zero, weilPairingVal_zero_left, mul_one]
      | Point.some x₂ y₂ hP₂, h₂ =>
        obtain ⟨xn, yn, hns, hμP, μ, hμx, hμy, hμF⟩ :=
          exists_mulByNHom_alg_charZero W n hn
        rcases eq_or_ne (Point.some x₁ y₁ hP₁ + Point.some x₂ y₂ hP₂) 0 with hsum | hsum
        · rw [hsum, weilPairingVal_zero_left]
          have hnn : Point.some x₂ y₂ hP₂
              = Point.some x₁ (W.negY x₁ y₁) ((nonsingular_neg x₁ y₁).mpr hP₁) :=
            (neg_eq_of_add_eq_zero_right hsum).symm
          rw [hnn]
          exact (weilPairingVal_mul_neg W (isUnit_two_of_charZero (F := F))
            (four_ne_zero_of_charZero (F := F)) n hn
            (natCast_ne_zero_of_charZero (F := F) n) μ hμF hns hμP hμx hμy hP₁ h₁ hQ hQt).symm
        · obtain ⟨x₃, y₃, hP₃, h3⟩ := exists_some_of_ne_zero W _ hsum
          rw [h3]
          exact (weilPairingVal_mul W (isUnit_two_of_charZero (F := F))
            (four_ne_zero_of_charZero (F := F)) n hn
            (natCast_ne_zero_of_charZero (F := F) n) μ hμF hns hμP hμx hμy hP₁ hP₂ hP₃ h3
            h₁ h₂ hQ hQt).symm

/-- ★★★★★★★★★**交代性(標数 0・完全版)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem weilPairingVal_alt (n : ℕ) (hn : 1 ≤ n) (P : W.Point) (hP : n • P = 0) :
    weilPairingVal W n P P = 1 := by
  haveI := isDedekindDomain_coordinateRing (isUnit_two_of_charZero (F := F)) W
  match P, hP with
  | 0, _ => exact weilPairingVal_zero_left W n 0
  | Point.some x y hPx, hP =>
      obtain ⟨xn, yn, hns, hμP, μ, hμx, hμy, hμF⟩ := exists_mulByNHom_alg_charZero W n hn
      exact weilPairingVal_self W (isUnit_two_of_charZero (F := F))
        (four_ne_zero_of_charZero (F := F)) n hn
        (fun k hk1 _ => Nat.cast_ne_zero.2 (by omega)) μ hμF hns hμP hμx hμy hPx hP

end CharZero

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal_mul_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性——e_n(P,Q)·e_n(−P,Q) = 1)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_add_left.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_alt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
