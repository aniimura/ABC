import ABC3.Found.GaloisRep.OmegaOne

/-!
# Galois (G1) 第 11 ブロック —— ★★★**`omegaNum = 0` を `ψ` だけの恒等式に**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★帰納段の**形が確定した**

    omegaNum n = 0
      ⟺ psiComp n = ψ_n (a₁ (X ψ_n² − ψ_{n+1} ψ_{n-1}) + a₃ ψ_n²)

★`φ_n = X ψ_n² − ψ_{n+1} ψ_{n-1}` は**定義**(`rfl`)なので、
★★右辺は **`ψ` と係数だけ**になった。

## ★★左辺も `ψ` で書ける(第 7)

    psiComp n × ψ₂ = ψ(n-1)² ψ(n+2) − ψ(n-2) ψ(n+1)²

★★★したがって、両辺に `ψ₂` を掛ければ**完全に `ψ` の恒等式**である:

    ψ(n-1)² ψ(n+2) − ψ(n-2) ψ(n+1)²
      = ψ_n (a₁(X ψ_n² − ψ_{n+1}ψ_{n-1}) + a₃ ψ_n²) (a₁X + a₃)      （標数 2)

★これを `normEDS_even` / `normEDS_odd` で帰納する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaNum_eq_zero_iff` | ★★`omegaNum = 0` の同値変形 |
| `phi_def` | ★`φ_n` の定義(`rfl`) |
| `omegaNum_eq_zero_iff'` | ★★★**`ψ` だけの恒等式** |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★★**`omegaNum = 0` の同値変形**。 -/
theorem omegaNum_eq_zero_iff (n : ℤ) :
    omegaNum W n = 0
      ↔ psiComp W n = W.ψ n * (C (C W.a₁) * W.φ n + C (C W.a₃) * W.ψ n ^ 2) := by
  rw [omegaNum, sub_eq_zero]

/-- ★**`φ_n = X ψ_n² − ψ_{n+1} ψ_{n-1}`**(定義、`rfl`)。 -/
theorem phi_def (n : ℤ) :
    W.φ n = C X * W.ψ n ^ 2 - W.ψ (n + 1) * W.ψ (n - 1) := rfl

/-- ★★★**`omegaNum = 0` は `ψ` だけの恒等式である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★第 7(`psiComp × ψ₂`)と合わせると、両辺が完全に `ψ` で書ける
——`normEDS_even` / `normEDS_odd` で帰納できる形である。 -/
theorem omegaNum_eq_zero_iff' (n : ℤ) :
    omegaNum W n = 0
      ↔ psiComp W n = W.ψ n * (C (C W.a₁) * (C X * W.ψ n ^ 2 - W.ψ (n + 1) * W.ψ (n - 1))
          + C (C W.a₃) * W.ψ n ^ 2) := by
  rw [omegaNum_eq_zero_iff, phi_def]

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_eq_zero_iff'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——omegaNum = 0 を ψ だけの恒等式に)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
