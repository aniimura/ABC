import ABC3.Found.GaloisRep.TateLimit

/-!
# Galois (G2) 第 77 ブロック —— **★★★★★★★★G2 達成 `T_l E ≅ ℤ_l²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これで原文の `GL₂(ℤ_l)` の行き先が書けた

    T_l E = lim_n E[l^n]  ≃+  ℤ_l × ℤ_l

★★`Interface` 側の `tateModule` は**逆極限として定義**してあるので
(第 76 の `limTors` と定義的に等しい)、posit は同型 1 本だけである。

## ★★標数 0 を課した(逸脱の記録)

(G1) の個数勘定が `∀ k ≤ n, (k : K) ≠ 0` を要求するため、
`n = l^m` を `m → ∞` で走らせるには **標数 0** が要る。
★ABC は `ℚ` の代数閉包の上で使うので十分である(§9-398 の記録の延長)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `torsion_padic_addEquiv` | ★★★★★★★★**`T_l E ≃+ ℤ_l²`** |
| `tateModuleDataWitness` | ★`TateModuleData` の実装 |
| `TateModuleData.nonvacuous` | ★★★★★★★★**G2 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine ABC3.Interface.GaloisRep

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★★★**`T_l E ≃+ ℤ_l²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem torsion_padic_addEquiv [IsAlgClosed F] [CharZero F] (hΔ : W.Δ ≠ 0) (l : ℕ)
    [Fact l.Prime] :
    Nonempty (limTors W.toAffine.Point l ≃+ (ℤ_[l] × ℤ_[l])) := by
  refine addEquiv_limTors ?_ ?_ l
  · intro m hm
    exact (finite_torsion W m hm (Nat.cast_ne_zero.2 (by omega))).to_subtype
  · intro m hm
    exact torsion_card W hΔ m hm (fun k hk1 _ => Nat.cast_ne_zero.2 (by omega))

/-- ★★★★★★★★**`TateModuleData` の実装**。 -/
noncomputable def tateModuleDataWitness : TateModuleData where
  toTorsionStructureData := torsionStructureDataWitness
  freeRankTwo := by
    intro K _ _ _ _ W hell l hp
    haveI := hell
    haveI := hp
    have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
    exact torsion_padic_addEquiv W hΔ l

/-- ★★★★★★★★**`TateModuleData` は非空虚である**——G2 達成。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが Galois 表現論の 8 件のうち **G2** である。 -/
theorem TateModuleData.nonvacuous : Nonempty TateModuleData :=
  ⟨tateModuleDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def torsion_padic_addEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進 Tate 加群 T_l E ≅ Z_l²)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
