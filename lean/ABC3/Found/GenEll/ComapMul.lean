import ABC3.Found.GenEll.ComapAffine

/-!
# [GenEll] Proposition 1.4, (i) —— **`comap` は積を保つ**(一般)(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★公式を計算せずに済む道

`ComapAffine.lean` で**アフィンの場合**は取れた。一般の場合には
「`(I.comap f).ideal U = Ideal.map (f.appLE W U h) (I.ideal W)`」という
**公式**を経由すると `Γ(U.toScheme, ⊤) ≅ Γ(X, U)` の共役の帳簿が要る。

★★★**公式は要らない。** 必要なのは **2 つの事実**だけである:

1. **開埋め込みに沿った `comap` は積を保つ**(`comap_mul_of_isOpenImmersion`)
   ——`appIso` が同型なので `Ideal.comap_symm` + `Ideal.map_mul` で出る
2. **アフィン射に沿った `comap` は積を保つ**(`ComapAffine.lean`、取得済)

★★あとは `U.ι ≫ f = f.resLE W U h ≫ W.ι` に沿って `comap_comp` を往復するだけ——
**書き換えの連鎖で閉じる。**

## ★★被覆(ヘッダ)

各点 `x : X` について `f(x)` を含むアフィン `W ⊆ Y` を取り、
`f ⁻¹ᵁ W` の中に `x` を含むアフィン `U` を取る。
★これらは `X` を覆うので `ext_of_iSup_eq_top` が使える。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory CategoryTheory.Limits

variable {X Y : Scheme}

/-! ## ★★開埋め込みに沿った `comap` -/

/-- ★★**開埋め込みに沿った `comap` は積を保つ**。

★機構: `ideal_comap_of_isOpenImmersion` が `Ideal.comap (appIso).inv` を与え、
`appIso` は**同型**なので `Ideal.comap_symm` で `Ideal.map` に直せる。
そこで `Ideal.map_mul` が効く。 -/
theorem comap_mul_of_isOpenImmersion (I J : Y.IdealSheafData) (ι : X ⟶ Y)
    [IsOpenImmersion ι] :
    (I * J).comap ι = I.comap ι * J.comap ι := by
  refine Scheme.IdealSheafData.ext (funext fun U => ?_)
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion]
  show _ = (I.comap ι).ideal U * (J.comap ι).ideal U
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion,
    Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion]
  set e := (Scheme.Hom.appIso ι U.1).commRingCatIsoToRingEquiv with he
  have key : ∀ K : Ideal Γ(Y, ι ''ᵁ U.1),
      Ideal.comap ((Scheme.Hom.appIso ι U.1).inv).hom K = K.map (e : _ →+* _) :=
    fun K => Ideal.comap_symm e
  rw [key, key, key, Scheme.IdealSheafData.ideal_mul, Pi.mul_apply, Ideal.map_mul]

/-! ## ★★★局所での一致 -/

/-- ★★★**`f(U) ⊆ W` なるアフィン開集合の対の上では `comap` は積を保つ**。

★機構は `U.ι ≫ f = f.resLE W U h ≫ W.ι` に沿った `comap_comp` の往復である:

    `(I·J)^*|_U = ((I·J)^{W.ι})^{resLE} = (I^{W.ι} · J^{W.ι})^{resLE}`
                `= (I^{W.ι})^{resLE} · (J^{W.ι})^{resLE} = (I^*|_U) · (J^*|_U)`

★★2 段目は開埋め込み(`comap_mul_of_isOpenImmersion`)、
3 段目は**アフィン射**(`comap_mul_of_isAffine`——`U.toScheme`, `W.toScheme` はアフィン)。 -/
theorem comap_mul_comap_ι (I J : Y.IdealSheafData) (f : X ⟶ Y)
    (U : X.affineOpens) (W : Y.affineOpens) (h : U.1 ≤ f ⁻¹ᵁ W.1) :
    ((I * J).comap f).comap U.1.ι = (I.comap f * J.comap f).comap U.1.ι := by
  haveI : IsAffine U.1.toScheme := U.2
  haveI : IsAffine W.1.toScheme := W.2
  have hfac : U.1.ι ≫ f = f.resLE W.1 U.1 h ≫ W.1.ι :=
    (Scheme.Hom.resLE_comp_ι f h).symm
  have hL : ((I * J).comap f).comap U.1.ι
      = ((I * J).comap W.1.ι).comap (f.resLE W.1 U.1 h) := by
    rw [← Scheme.IdealSheafData.comap_comp, hfac, Scheme.IdealSheafData.comap_comp]
  have hI : (I.comap f).comap U.1.ι
      = (I.comap W.1.ι).comap (f.resLE W.1 U.1 h) := by
    rw [← Scheme.IdealSheafData.comap_comp, hfac, Scheme.IdealSheafData.comap_comp]
  have hJ : (J.comap f).comap U.1.ι
      = (J.comap W.1.ι).comap (f.resLE W.1 U.1 h) := by
    rw [← Scheme.IdealSheafData.comap_comp, hfac, Scheme.IdealSheafData.comap_comp]
  rw [hL, comap_mul_of_isOpenImmersion, comap_mul_of_isAffine,
    comap_mul_of_isOpenImmersion, hI, hJ]

/-! ## ★★開埋め込みへの制限は情報を落とさない -/

/-- ★★**アフィン開集合への制限が一致すれば、その上での `ideal` も一致する**。

★★★**器具の要点**: 同型の逆を打ち消すのに `RingEquiv` を経由しない。
`Iso.hom_inv_id` で**圏の言葉のまま**打ち消す。
★`commRingCatIsoToRingEquiv` + `Ideal.comap_map_of_bijective` を使うと
`whnf` が 200 万ヒートビートでも終わらない(2026-08-17 実測)。 -/
theorem ideal_eq_of_comap_ι_eq {A B : X.IdealSheafData} (U : X.affineOpens)
    (h : A.comap U.1.ι = B.comap U.1.ι) : A.ideal U = B.ideal U := by
  haveI : IsAffine U.1.toScheme := U.2
  have hUeq : (⟨U.1.ι ''ᵁ (⊤ : U.1.toScheme.Opens),
      (isAffineOpen_top U.1.toScheme).image_of_isOpenImmersion U.1.ι⟩ : X.affineOpens) = U := by
    apply Subtype.ext; simp
  have h2 := congrArg (fun K => K.ideal ⟨⊤, isAffineOpen_top U.1.toScheme⟩) h
  simp only [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion] at h2
  have hcomp : ∀ K : Ideal Γ(X, U.1.ι ''ᵁ (⊤ : U.1.toScheme.Opens)),
      Ideal.comap ((Scheme.Hom.appIso U.1.ι (⊤ : U.1.toScheme.Opens)).hom).hom
        (Ideal.comap ((Scheme.Hom.appIso U.1.ι (⊤ : U.1.toScheme.Opens)).inv).hom K) = K := by
    intro K
    rw [Ideal.comap_comap, ← CommRingCat.hom_comp, Iso.hom_inv_id]
    simp
  have hinj := congrArg
    (Ideal.comap ((Scheme.Hom.appIso U.1.ι (⊤ : U.1.toScheme.Opens)).hom).hom) h2
  rw [hcomp, hcomp] at hinj
  rw [← hUeq]
  exact hinj

/-! ## ★★被覆 -/

/-- ★**各点で「`f` の像がアフィンに収まるアフィン近傍」が取れる**。

★`f(x)` を含むアフィン `W ⊆ Y` を取り、`f ⁻¹ᵁ W` の中に
`x` を含むアフィン `U` を取るだけ(`isBasis_affineOpens` を 2 回)。 -/
theorem exists_affinePair (f : X ⟶ Y) (x : X) :
    ∃ (W : Y.affineOpens) (U : X.affineOpens), x ∈ U.1 ∧ U.1 ≤ f ⁻¹ᵁ W.1 := by
  obtain ⟨_, ⟨W, hW, rfl⟩, hfx, -⟩ :=
    Y.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ (f.base x)) isOpen_univ
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, hUsub⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open
      (show x ∈ (f ⁻¹ᵁ W : X.Opens) from hfx) (f ⁻¹ᵁ W).2
  exact ⟨⟨W, hW⟩, ⟨U, hU⟩, hxU, hUsub⟩

/-! ## ★★★`comap` は積を保つ(一般) -/

/-- ★★★**`Scheme.IdealSheafData.comap` は積を保つ**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**これが `Proposition 1.4, (i)` に残っていた 1 本である。**
`HeightConstruction.lean` の `PullbackMul` がこれで**無条件に**満たされる。

★★機構は 3 段:
1. 各点で `f(U) ⊆ W` なるアフィンの対を取る(`exists_affinePair`)
2. その上では `comap` が積を保つ(`comap_mul_comap_ι`)
3. 被覆から貼る(`ext_of_iSup_eq_top`)

★★★**mathlib に無かったもの**(2026-08-17 実測。`comap_sup` はあるが `comap_mul` は無い)。 -/
theorem comap_mul (I J : Y.IdealSheafData) (f : X ⟶ Y) :
    (I * J).comap f = I.comap f * J.comap f := by
  choose W U hxU hUW using exists_affinePair f
  refine Scheme.IdealSheafData.ext_of_iSup_eq_top U ?_ (fun x => ?_)
  · refine top_le_iff.1 fun x _ => ?_
    exact TopologicalSpace.Opens.mem_iSup.2 ⟨x, hxU x⟩
  · exact ideal_eq_of_comap_ι_eq (U x) (comap_mul_comap_ι I J f (U x) (W x) (hUW x))

/-! ## ★出典の紐付け(`.src`) -/

def comap_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(引き戻しが積を保つこと)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
