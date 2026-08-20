import ABC3.Found.GaloisRep.WeilGalPoint

/-!
# Galois (G5) 第 195 ブロック —— **★★★★★★★★第 2 変数の双線型性と反対称性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`det ρ` の計算に要る道具

`det_galRep_eq_cyclotomic` は基底 `P, Q` を取って

    e_n(aP + cQ, bP + dQ) = e_n(P,Q)^{ad − bc}

を出す。★これには**第 2 変数の双線型性**と**反対称性**、そして `k` 倍の公式が要る。

### ★★★★★★第 2 変数は平行移動の合成から出る

`e_n(P,Q) = τ_Q(g)/g` で、第 189 の `isTranslate_trans` から
**`τ_{Q₁+Q₂} = τ_{Q₁} ∘ τ_{Q₂}`**(第 190 の一意性)。★したがって

    τ_{Q₁+Q₂}(g)/g = τ_{Q₁}(τ_{Q₂}(g))/g = τ_{Q₁}(c₂·g)/g = c₂·(τ_{Q₁}(g)/g) = c₁c₂

★★`τ_{Q₁}` が定数 `c₂` を固定することを使う(`F`-代数同型だから)。

### ★★★値そのものが `τ_Q(g)/g` である

`Q = O` の場合も込めて **`algebraMap(e_n(P,Q)) = τ_Q(g)/g`** を出しておくと、
上の計算がそのまま値の等式になる(`τ_O = id` なので `Q = O` では両辺 1)。

### ★★★★★★反対称性は 3 行

    1 = e(P+Q, P+Q) = e(P,P)·e(P,Q)·e(Q,P)·e(Q,Q) = e(P,Q)·e(Q,P)

★交代性(第 190)と両変数の双線型性から。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `translateAut_fix_mu` | `n·Q = O` なら `τ_Q` は `μ` の像を固定する(`Q = O` 込み) |
| `translateAut_div_add` | ★★★★★★**`τ_{Q₁+Q₂}(g)/g = (τ_{Q₁}g/g)(τ_{Q₂}g/g)`** |
| `algebraMap_weilPairingVal` | ★★★★★★★**値は `τ_Q(g)/g` そのもの** |
| `weilPairingVal_add_right` | ★★★★★★★★**第 2 変数の双線型性** |
| `weilPairingVal_antisymm` | ★★★★★★**反対称性** |
| `weilPairingVal_nsmul_left` / `_right` | ★★★★★**`e_n(kP,Q) = e_n(P,Q)^k`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★平行移動の側の計算 -/

section Translate

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-- `translateAut` を `τ` に取った witness で値が決まる。 -/
theorem weilPairingVal_of_translateAut (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {x y : F} (hP : W.Nonsingular x y) (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {g : W.FunctionField} (hg : g ^ n = μ fP)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) {c : F}
    (hd : translateAut W h4 (Point.some x₀ y₀ hQ) g / g = algebraMap F W.FunctionField c) :
    weilPairingVal W n (Point.some x y hP) (Point.some x₀ y₀ hQ) = c := by
  have hτ := isTranslate_translateAut W h4 (Point.some x₀ y₀ hQ)
  exact weilPairingVal_eq W h2 hn
    ⟨x, y, hP, x₀, y₀, hQ, fP, μ, g, translateAut W h4 (Point.some x₀ y₀ hQ), xn, yn, hns,
      rfl, rfl, hfP, hinj, hμF, hμx, hμy, hμP, hg, hτ.1, hτ.2, hd⟩

/-- ★★★`n·Q = O` なら `τ_Q` は `μ` の像を固定する(`Q = O` 込み)。 -/
theorem translateAut_fix_mu (h4 : (4 : F) ≠ 0) (n : ℕ)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    {Q : W.Point} (hQt : n • Q = 0) (r : W.CoordinateRing) :
    translateAut W h4 Q (μ r) = μ r := by
  have hT := isTranslate_translateAut W h4 Q
  match Q, hQt, hT with
  | 0, _, hT => rw [show translateAut W h4 (0 : W.Point) = AlgEquiv.refl from hT]; rfl
  | Point.some x₀ y₀ hQ, hQt, hT =>
      have hcomp := aut_comp_mulByN W (translateAut W h4 (Point.some x₀ y₀ hQ)) hQ hT.1 hT.2
        n hQt μ hμF hns hμx hμy hμP
      exact congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun)
        hcomp) r

/-- ★★★★★★**`τ_{Q₁+Q₂}(g)/g = (τ_{Q₁}(g)/g)·(τ_{Q₂}(g)/g)`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 190 の `translateAut_trans`(平行移動の合成)そのもの。 -/
theorem translateAut_div_add (h4 : (4 : F) ≠ 0) {g : W.FunctionField} (hg0 : g ≠ 0)
    {Q₁ Q₂ : W.Point} {c₁ c₂ : F}
    (hd₁ : translateAut W h4 Q₁ g / g = algebraMap F W.FunctionField c₁)
    (hd₂ : translateAut W h4 Q₂ g / g = algebraMap F W.FunctionField c₂) :
    translateAut W h4 (Q₁ + Q₂) g / g = algebraMap F W.FunctionField (c₁ * c₂) := by
  have h2g : translateAut W h4 Q₂ g = algebraMap F W.FunctionField c₂ * g :=
    (div_eq_iff hg0).1 hd₂
  have h1g : translateAut W h4 Q₁ g = algebraMap F W.FunctionField c₁ * g :=
    (div_eq_iff hg0).1 hd₁
  have hcomp : translateAut W h4 (Q₁ + Q₂) g
      = translateAut W h4 Q₁ (translateAut W h4 Q₂ g) := by
    rw [← translateAut_trans W h4 Q₁ Q₂]
    rfl
  rw [hcomp, h2g, map_mul, (translateAut W h4 Q₁).commutes, h1g, map_mul]
  field_simp

/-- ★★★★★★★**値は `τ_Q(g)/g` そのものである**(`Q = O` 込み)。 -/
theorem algebraMap_weilPairingVal (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {x y : F} (hP : W.Nonsingular x y) (fP : W.CoordinateRing) (hfP0 : fP ≠ 0)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {g : W.FunctionField} (hg : g ^ n = μ fP)
    {Q : W.Point} (hQt : n • Q = 0) :
    algebraMap F W.FunctionField (weilPairingVal W n (Point.some x y hP) Q)
      = translateAut W h4 Q g / g := by
  have hfne : μ fP ≠ 0 := fun h0 => hfP0 (hinj (by rw [h0, map_zero]))
  have hg0 : g ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hg
    exact hfne hg.symm
  have hfix := translateAut_fix_mu W h4 n μ hμF hns hμx hμy hμP hQt fP
  obtain ⟨c, hc0, hd⟩ := exists_const_aut_div W h2 hn hg hfne hfix
  rw [hd]
  congr 1
  match Q, hQt, hd with
  | 0, _, hd =>
      rw [weilPairingVal_zero_right]
      rw [show translateAut W h4 (0 : W.Point) = AlgEquiv.refl from
        isTranslate_translateAut W h4 0] at hd
      have hone : (1 : W.FunctionField) = algebraMap F W.FunctionField c := by
        rw [← hd]; exact (div_self hg0).symm
      exact ((algebraMap F W.FunctionField).injective (by rw [← hone, map_one])).symm
  | Point.some x₀ y₀ hQ, hQt, hd =>
      exact weilPairingVal_of_translateAut W h2 h4 n hn hinj hμF hns hμP hμx hμy hP fP hfP hg
        hQ hd

end Translate

/-! ## ★★★★★★★★第 2 変数の双線型性 -/

section CharZero

variable {F : Type} [Field F] [DecidableEq F] [CharZero F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★★★**双線型性(第 2 変数、標数 0・完全版)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★平行移動の合成 `τ_{Q₁+Q₂} = τ_{Q₁} ∘ τ_{Q₂}` から出る。 -/
theorem weilPairingVal_add_right (n : ℕ) (hn : 1 ≤ n) (P Q₁ Q₂ : W.Point)
    (hPt : n • P = 0) (h₁ : n • Q₁ = 0) (h₂ : n • Q₂ = 0) :
    weilPairingVal W n P (Q₁ + Q₂)
      = weilPairingVal W n P Q₁ * weilPairingVal W n P Q₂ := by
  have h2' : IsUnit (2 : F) := Ne.isUnit (two_ne_zero (α := F))
  have h4' : (4 : F) ≠ 0 := by
    have h : ((4 : ℕ) : F) ≠ 0 := Nat.cast_ne_zero.2 (by norm_num)
    simp only [Nat.cast_ofNat] at h
    exact h
  have hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0 :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  haveI := isDedekindDomain_coordinateRing h2' W
  match P, hPt with
  | 0, _ => rw [weilPairingVal_zero_left, weilPairingVal_zero_left, weilPairingVal_zero_left,
      one_mul]
  | Point.some x y hP, hPt =>
      obtain ⟨xn, yn, hns, hμP, μ, hμx, hμy, hμF⟩ := exists_mulByNHom_alg_charZero W n hn
      have hinj := mulByN_injective W h2' n hn (hchar n hn le_rfl) μ hμF hns hμx hμy hμP
      obtain ⟨fP, hfP0, hfP⟩ := xyIdeal_pow_isPrincipal_integral hP n hPt
      obtain ⟨g, hg⟩ := exists_nthRoot_mulByN W h2' hP n hn hchar hPt μ hinj hμF hns hμP
        hμx hμy fP hfP
      have hfne : μ fP ≠ 0 := fun h0 => hfP0 (hinj (by rw [h0, map_zero]))
      have hg0 : g ≠ 0 := by
        intro h0
        rw [h0, zero_pow (by omega : n ≠ 0)] at hg
        exact hfne hg.symm
      have hkey := fun {Q : W.Point} (hQ : n • Q = 0) =>
        algebraMap_weilPairingVal W h2' h4' n hn hinj hμF hns hμP hμx hμy hP fP hfP0 hfP hg hQ
      refine (algebraMap F W.FunctionField).injective ?_
      rw [hkey (show n • (Q₁ + Q₂) = 0 by rw [nsmul_add, h₁, h₂, add_zero]),
        translateAut_div_add W h4' hg0 (hkey h₁).symm (hkey h₂).symm]

/-- ★★★★★★**反対称性** `e_n(P,Q)·e_n(Q,P) = 1`。 -/
theorem weilPairingVal_antisymm (n : ℕ) (hn : 1 ≤ n) (P Q : W.Point)
    (hP : n • P = 0) (hQ : n • Q = 0) :
    weilPairingVal W n P Q * weilPairingVal W n Q P = 1 := by
  have hPQ : n • (P + Q) = 0 := by rw [nsmul_add, hP, hQ, add_zero]
  have h1 := weilPairingVal_alt W n hn (P + Q) hPQ
  rw [weilPairingVal_add_left W n hn P Q (P + Q) hP hQ hPQ,
    weilPairingVal_add_right W n hn P P Q hP hP hQ,
    weilPairingVal_add_right W n hn Q P Q hQ hP hQ,
    weilPairingVal_alt W n hn P hP, weilPairingVal_alt W n hn Q hQ, one_mul, mul_one] at h1
  exact h1

/-- ★★★★★**`e_n(k·P, Q) = e_n(P,Q)^k`**。 -/
theorem weilPairingVal_nsmul_left (n : ℕ) (hn : 1 ≤ n) (P Q : W.Point)
    (hP : n • P = 0) (hQ : n • Q = 0) (k : ℕ) :
    weilPairingVal W n (k • P) Q = (weilPairingVal W n P Q) ^ k := by
  induction k with
  | zero => rw [zero_smul, pow_zero, weilPairingVal_zero_left]
  | succ k ih =>
      have hk : n • (k • P) = 0 := by rw [smul_comm, hP, smul_zero]
      rw [succ_nsmul, weilPairingVal_add_left W n hn (k • P) P Q hk hP hQ, ih, pow_succ]

/-- ★★★★★**`e_n(P, k·Q) = e_n(P,Q)^k`**。 -/
theorem weilPairingVal_nsmul_right (n : ℕ) (hn : 1 ≤ n) (P Q : W.Point)
    (hP : n • P = 0) (hQ : n • Q = 0) (k : ℕ) :
    weilPairingVal W n P (k • Q) = (weilPairingVal W n P Q) ^ k := by
  induction k with
  | zero => rw [zero_smul, pow_zero, weilPairingVal_zero_right]
  | succ k ih =>
      have hk : n • (k • Q) = 0 := by rw [smul_comm, hQ, smul_zero]
      rw [succ_nsmul, weilPairingVal_add_right W n hn P (k • Q) Q hP hk hQ, ih, pow_succ]

end CharZero

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal_add_right.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性——第 2 変数)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_antisymm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の反対称性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
