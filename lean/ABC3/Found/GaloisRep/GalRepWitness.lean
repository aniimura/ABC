import ABC3.Found.GaloisRep.GalRep

/-!
# Galois (G3) 第 92 ブロック —— **★★★★★★★★G3 達成**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★原文が名指しする表現が obligation として閉じた

    ρ_{E,l} : Gal(L/K) → GL₂(ℤ_l)

★これは `T_l E` への**本物の Galois 作用**の行列表示である
(`rep_action`——空虚な posit ではない)。

## ★★★`det = 円分指標` は (G5) へ移した(§9-407)

原文がこの事実を使うのは**像の主張の側**である——
`det` が全射だからこそ「像が `GL₂` 全体」まで言える。
★★内容は **Weil 対**(mathlib に 0 件)であり、(G3) の構成とは独立の定理である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `nonempty_tate_basis` | ★`T_l E` の基底の存在((G2) から) |
| `repChoice` | ★基底を選んで得る表現 |
| `galoisRepDataWitness` | ★`GaloisRepData` の実装 |
| `GaloisRepData.nonvacuous` | ★★★★★★★★**G3 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★`T_l E` の基底が存在する((G2) の帰結)。 -/
theorem nonempty_tate_basis [IsAlgClosed L] [CharZero K] (W : WeierstrassCurve K)
    (hell : W.IsElliptic) (l : ℕ) [Fact l.Prime] :
    Nonempty (tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) := by
  haveI : CharZero L := charZero_of_injective_algebraMap (algebraMap K L).injective
  haveI hE : (W.baseChange L).IsElliptic := by
    haveI := hell
    unfold WeierstrassCurve.baseChange
    infer_instance
  have hΔ : (W.baseChange L).Δ ≠ 0 := (W.baseChange L).isUnit_Δ.ne_zero
  obtain ⟨e0⟩ := torsion_padic_addEquiv (W.baseChange L) hΔ l
  exact ⟨e0.trans (prodEquivPiTwo l)⟩

open scoped Classical in
/-- ★基底を 1 つ選んで得る表現(基底が無いところでは自明表現)。 -/
noncomputable def repChoice (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] :
    (L ≃ₐ[K] L) →* Matrix.GeneralLinearGroup (Fin 2) ℤ_[l] :=
  if h : Nonempty (tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
  then galRep W l (Classical.choice h)
  else 1

open scoped Classical in
theorem repChoice_eq (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (h : Nonempty (tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))) :
    repChoice W l = galRep W l (Classical.choice h) := by
  rw [repChoice, dif_pos h]

/-- ★★★★★★★★**`GaloisRepData` の実装**。 -/
noncomputable def galoisRepDataWitness : GaloisRepData where
  toTateModuleData := tateModuleDataWitness
  rep := by
    intro K L _ _ _ _ _ W l hp σ
    haveI := hp
    exact repChoice W l σ
  rep_mul := by
    intro K L _ _ _ _ _ W l hp σ τ
    haveI := hp
    exact map_mul (repChoice W l) σ τ
  rep_action := by
    intro K L _ _ _ _ _ _ _ W hell l hp
    haveI := hp
    have hne : Nonempty (tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :=
      nonempty_tate_basis W hell l
    refine ⟨Classical.choice hne, fun σ x => ?_⟩
    have hrep : repChoice W l = galRep W l (Classical.choice hne) := repChoice_eq W l hne
    rw [hrep, galRep_coe]
    exact galMat_apply W l _ σ x

/-- ★★★★★★★★**`GaloisRepData` は非空虚である**——G3 達成。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが Galois 表現論の 8 件のうち **G3** である。 -/
theorem GaloisRepData.nonvacuous : Nonempty GaloisRepData :=
  ⟨galoisRepDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def galoisRepDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進 Galois 表現の構成)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
