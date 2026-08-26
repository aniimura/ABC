import ABC3.Found.GaloisRep.WeierstrassAddition

/-!
# Galois (G6) 第 249 ブロック —— **★★★★★★★★解析側の共線性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★葉 (c) 段 3 の入口 —— `℘` の共線性を Tate 座標に移す

第 248 の `weierstrass_collinear` を、第 218・220 の対応

    ℘(z)  = (2πi)²(X + 1/12)          (`weierstrassP_eq_tateXfun`)
    ℘'(z) = (2πi)³(2Y + X)            (`derivWeierstrassP_eq_tateYfun`)

で書き換えると、**`1/12` の項も `XᵢXⱼ` の項もちょうど消えて**

    X₁(Y₂ − Y₃) + X₂(Y₃ − Y₁) + X₃(Y₁ − Y₂) = 0

という**同じ形の行列式**が残る(`collinear_tatefun`)。★係数は `2(2πi)⁵` だけである。

## ★★★★★対称な 3 変数 `(u, v, w)` へ

`z₃ := τ − z₁ − z₂` と取ると `z₁ + z₂ + z₃ = τ` なので、
`u := e^{2πi z₁}`、`v := e^{2πi z₂}`、`w := e^{2πi z₃}` について

    q = u v w

となる。しかも各点の「相方」`q/uᵢ` が

    q/u = v w,   q/v = u w,   q/w = u v

と**すべて多項式**になる。★★★これが段 3 の普遍環を `ℤ[U,V,W]`(`q = UVW`)に
取れる理由である——第 223 の `ℤ[A,W]`(`q = AW`)の素直な 3 変数化になる。

## ★上半平面から自由になる

第 237 の `exists_uh_exp_eq`(`‖t‖ < 1`、`t ≠ 0` なら `t = e^{2πi z}` となる `z ∈ ℍ` がある)
で `u, v, w` を任意の「単位開円板の 0 でない点」に取り替えられる。
★`τ := z₁ + z₂ + z₃` と**定義してしまう**ので、`im` の条件は `z₃ ∈ ℍ` から自動的に出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `thirdUH` | ★第 3 の点 `τ − z₁ − z₂` |
| `collinear_tatefun` | ★★★★★★★`tateXfun`/`tateYfun` の形の共線性 |
| `collinear_analytic_uvw` | ★★★★★★★★**`(u,v,w)` の形の共線性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real Filter Topology PeriodPair

/-! ## ★第 3 の点 -/

/-- ★`im z₁ + im z₂ < im τ` なら `τ − z₁ − z₂` も上半平面の点。 -/
noncomputable def thirdUH (z₁ z₂ τ : UpperHalfPlane)
    (h : (z₁ : ℂ).im + (z₂ : ℂ).im < (τ : ℂ).im) : UpperHalfPlane :=
  ⟨(τ : ℂ) - (z₁ : ℂ) - (z₂ : ℂ), by
    simp only [Complex.sub_im]
    linarith⟩

@[simp] theorem thirdUH_coe (z₁ z₂ τ : UpperHalfPlane)
    (h : (z₁ : ℂ).im + (z₂ : ℂ).im < (τ : ℂ).im) :
    ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ) = (τ : ℂ) - (z₁ : ℂ) - (z₂ : ℂ) := rfl

/-! ## ★★★★★★★`tateXfun`/`tateYfun` の形 -/

/-- ★★★★★★★**解析側の共線性(`tateXfun`/`tateYfun` の形)**。

`℘ = (2πi)²(X + 1/12)`、`℘' = (2πi)³(2Y + X)` を第 248 に代入すると、
`1/12` の項も `XᵢXⱼ` の項も消えて **`2(2πi)⁵` 倍の同じ行列式**が残る。 -/
theorem collinear_tatefun (z₁ z₂ τ : UpperHalfPlane)
    (h : (z₁ : ℂ).im + (z₂ : ℂ).im < (τ : ℂ).im) :
    tateXfun z₁ τ * (tateYfun z₂ τ - tateYfun (thirdUH z₁ z₂ τ h) τ)
      + tateXfun z₂ τ * (tateYfun (thirdUH z₁ z₂ τ h) τ - tateYfun z₁ τ)
      + tateXfun (thirdUH z₁ z₂ τ h) τ * (tateYfun z₁ τ - tateYfun z₂ τ) = 0 := by
  have hp1 : (0 : ℝ) < (z₁ : ℂ).im := z₁.2
  have hp2 : (0 : ℝ) < (z₂ : ℂ).im := z₂.2
  have him1 : (z₁ : ℂ).im < (τ : ℂ).im := by linarith
  have him2 : (z₂ : ℂ).im < (τ : ℂ).im := by linarith
  have him3 : ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ).im < (τ : ℂ).im := by
    rw [thirdUH_coe]
    simp only [Complex.sub_im]
    linarith
  have hs : (z₁ : ℂ) + (z₂ : ℂ) + ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ)
      ∈ (tauPair τ).lattice := by
    rw [thirdUH_coe, show (z₁ : ℂ) + (z₂ : ℂ) + ((τ : ℂ) - (z₁ : ℂ) - (z₂ : ℂ)) = (τ : ℂ) by ring]
    exact (tauPair τ).ω₁_mem_lattice
  have hkey := weierstrass_collinear (tauPair τ) (z₁ : ℂ) (z₂ : ℂ)
    ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ) hs
    (notMem_lattice_of_im_lt z₁ τ him1) (notMem_lattice_of_im_lt z₂ τ him2)
    (notMem_lattice_of_im_lt (thirdUH z₁ z₂ τ h) τ him3)
  rw [weierstrassP_eq_tateXfun z₁ τ him1, weierstrassP_eq_tateXfun z₂ τ him2,
    weierstrassP_eq_tateXfun (thirdUH z₁ z₂ τ h) τ him3,
    derivWeierstrassP_eq_tateYfun z₁ τ him1, derivWeierstrassP_eq_tateYfun z₂ τ him2,
    derivWeierstrassP_eq_tateYfun (thirdUH z₁ z₂ τ h) τ him3] at hkey
  have hne : (2 : ℂ) * (2 * ↑π * I) ^ 5 ≠ 0 :=
    mul_ne_zero two_ne_zero (pow_ne_zero _ two_pi_I_ne_zero)
  refine (mul_eq_zero.1 ?_).resolve_left hne
  linear_combination hkey

/-! ## ★★★★★★★★`(u,v,w)` の形 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**解析側の共線性(`(u,v,w)` の形)**——3 点は `u, v, w`、`q = u v w`。

各点の相方 `q/uᵢ` が `v w`, `u w`, `u v` と**すべて多項式**になるのが要点である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collinear_analytic_uvw (u v w : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0) (hw0 : w ≠ 0)
    (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    tateXanalytic u (v * w) * (tateYanalytic v (u * w) - tateYanalytic w (u * v))
      + tateXanalytic v (u * w) * (tateYanalytic w (u * v) - tateYanalytic u (v * w))
      + tateXanalytic w (u * v) * (tateYanalytic u (v * w) - tateYanalytic v (u * w)) = 0 := by
  obtain ⟨z₁, hz₁⟩ := exists_uh_exp_eq u hu0 hu
  obtain ⟨z₂, hz₂⟩ := exists_uh_exp_eq v hv0 hv
  obtain ⟨z₃, hz₃⟩ := exists_uh_exp_eq w hw0 hw
  have hp1 : (0 : ℝ) < (z₁ : ℂ).im := z₁.2
  have hp2 : (0 : ℝ) < (z₂ : ℂ).im := z₂.2
  have hp3 : (0 : ℝ) < (z₃ : ℂ).im := z₃.2
  let τ : UpperHalfPlane := ⟨(z₁ : ℂ) + (z₂ : ℂ) + (z₃ : ℂ), by
    simp only [Complex.add_im]; linarith⟩
  have hτ : (τ : ℂ) = (z₁ : ℂ) + (z₂ : ℂ) + (z₃ : ℂ) := rfl
  have h : (z₁ : ℂ).im + (z₂ : ℂ).im < (τ : ℂ).im := by
    rw [hτ]; simp only [Complex.add_im]; linarith
  have him1 : (z₁ : ℂ).im < (τ : ℂ).im := by rw [hτ]; simp only [Complex.add_im]; linarith
  have him2 : (z₂ : ℂ).im < (τ : ℂ).im := by rw [hτ]; simp only [Complex.add_im]; linarith
  have him3 : ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ).im < (τ : ℂ).im := by
    rw [thirdUH_coe]; simp only [Complex.sub_im]; linarith
  have hc3 : ((thirdUH z₁ z₂ τ h : UpperHalfPlane) : ℂ) = (z₃ : ℂ) := by
    rw [thirdUH_coe, hτ]; ring
  have hkey := collinear_tatefun z₁ z₂ τ h
  rw [tateXfun_eq_analytic z₁ τ him1, tateXfun_eq_analytic z₂ τ him2,
    tateXfun_eq_analytic (thirdUH z₁ z₂ τ h) τ him3,
    tateYfun_eq_analytic z₁ τ him1, tateYfun_eq_analytic z₂ τ him2,
    tateYfun_eq_analytic (thirdUH z₁ z₂ τ h) τ him3] at hkey
  have e1 : Complex.exp (2 * ↑π * I * ((subUH z₁ τ him1 : UpperHalfPlane) : ℂ)) = v * w := by
    rw [subUH_coe, hτ, ← hz₂, ← hz₃, ← Complex.exp_add]
    ring_nf
  have e2 : Complex.exp (2 * ↑π * I * ((subUH z₂ τ him2 : UpperHalfPlane) : ℂ)) = u * w := by
    rw [subUH_coe, hτ, ← hz₁, ← hz₃, ← Complex.exp_add]
    ring_nf
  have e3 : Complex.exp (2 * ↑π * I *
      ((subUH (thirdUH z₁ z₂ τ h) τ him3 : UpperHalfPlane) : ℂ)) = u * v := by
    rw [subUH_coe, hc3, hτ, ← hz₁, ← hz₂, ← Complex.exp_add]
    ring_nf
  rw [hz₁, hz₂, hc3, hz₃, e1, e2, e3] at hkey
  linear_combination hkey

/-! ## ★出典の紐付け(`.src`) -/

def collinear_tatefun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の共線性、tateXfun の形)",
    sectionId := "genell-def-3-3" }

def collinear_analytic_uvw.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の共線性、(u,v,w) の形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
