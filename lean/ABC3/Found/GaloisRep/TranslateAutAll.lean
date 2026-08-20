import ABC3.Found.GaloisRep.NotTwoTorsionPoint

/-!
# Galois (G5) 第 124 ブロック —— **★★★★★★葉 1(平行移動)が完全に閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★2 等分点でも自己同型になる

第 120 で 2 等分点でない `Q` の自己同型を得た。★2 等分点の場合は
第 122 の分解 `Q = Q₁ + Q₂`(どちらも 2 等分点でない)を使い、

    τ_Q = τ_{Q₁} ∘ τ_{Q₂}

を**自己同型の合成**として作る。★★合成則(第 121)が座標の値を保証する。

## ★★★★★これで (G5) の葉 1 は終わりである

    平行移動 τ_Q : F(W) ≃ₐ[F] F(W)  がすべての Q に対して存在する

★第 114-124 の **11 ブロック**で閉じた(当初の見積もり 20-60)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_translateAut_all` | ★★★★★★**すべての `Q` で自己同型** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★★★★**平行移動はつねに関数体の `F` 自己同型を誘導する**
(代数閉・標数 ≠ 2)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★2 等分点でない場合は第 120。★★2 等分点の場合は第 122 の分解を使って
自己同型の合成として作り、第 121 の合成則で座標の値を出す。 -/
theorem exists_translateAut_all [IsAlgClosed F] [Infinite F] [DecidableEq F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (h4 : (4 : F) ≠ 0)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField,
      τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀ := by
  by_cases h2 : W.negY x₀ y₀ = y₀
  · obtain ⟨x₁, y₁, x₂, y₂, h₁, h₂, k₁, k₂, hsum⟩ :=
      exists_decomp_of_twoTorsion hQ h2 (exists_notTwoTorsion_point W h4)
    have hmap := congrArg (toFF W) hsum
    rw [map_add, toFF_some, toFF_some, toFF_some] at hmap
    obtain ⟨hx, hy⟩ := translateAlgHom_comp_general W h₁ h₂ hQ k₁ hmap
    refine ⟨(translateAlgEquiv W h₂ k₂).trans (translateAlgEquiv W h₁ k₁), ?_, ?_⟩
    · rw [AlgEquiv.trans_apply, translateAlgEquiv_coordX]
      exact hx
    · rw [AlgEquiv.trans_apply, translateAlgEquiv_coordY]
      exact hy
  · exact exists_translateAut_of_not_twoTorsion W hQ h2

/-! ## ★出典の紐付け(`.src`) -/

def exists_translateAut_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動がつねに関数体の自己同型であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
