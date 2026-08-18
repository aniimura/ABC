import ABC3.Found.Arakelov.PicGammaUnit

/-!
# Arakelov (B1) 第 69 ブロック —— **係数環を上げるテンソルの射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★比較射の前段

第 66 ブロックで比較射の形が

    F.val(⊤) ⊗_R G.val(⊤)  →  (F.val ⊗ G.val)(⊤)  →  (層化 ..)(⊤)
    ───── 前段(本ブロック) ──   ── 後段(第 68) ──

と分かった。★**前段は係数環を `R` から `𝒪(⊤)` へ上げる射**である。

## ★★機構

`R → S` と `S` 加群 `M, N` に対し

    M ⊗[R] N  →  M ⊗[S] N,   m ⊗ₜ n ↦ m ⊗ₜ n

は `TensorProduct.lift` で作れる(mathlib に無いので書く、一発で通った)。

★★★第 67 ブロックで `Γ(Spec R, ⊤) ≅ R` を確認してあるので、
**この射は同型になる**(係数環が同型だから)。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tensorBaseUp` | ★★★**`M ⊗[R] N → M ⊗[S] N`**(汎用) |
| `tensorBaseUp_tmul` | ★純テンソルでの値 |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory

section

variable {Rr Ss : Type u} [CommRing Rr] [CommRing Ss] [Algebra Rr Ss]
  (M N : Type u) [AddCommGroup M] [Module Rr M] [Module Ss M] [IsScalarTower Rr Ss M]
  [AddCommGroup N] [Module Rr N] [Module Ss N] [IsScalarTower Rr Ss N]

/-- ★★★**係数環を上げるテンソルの射** `M ⊗[R] N → M ⊗[S] N`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが `Γ(F) ⊗ Γ(G) → Γ(F ⊗ G)` の前段である。 -/
noncomputable def tensorBaseUp :
    TensorProduct Rr M N →ₗ[Rr] TensorProduct Ss M N :=
  TensorProduct.lift
    { toFun := fun m => (TensorProduct.mk Ss M N m).restrictScalars Rr
      map_add' := by intro m m'; ext n; simp
      map_smul' := by intro r m; ext n; simp }

/-- ★**純テンソルでの値**。 -/
@[simp] theorem tensorBaseUp_tmul (m : M) (n : N) :
    tensorBaseUp (Rr := Rr) (Ss := Ss) M N (m ⊗ₜ[Rr] n) = m ⊗ₜ[Ss] n := rfl

end

/-! ## ★出典の紐付け(`.src`) -/

def tensorBaseUp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数環を上げるテンソルの射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
