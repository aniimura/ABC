/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DescentUnique

/-!
# [GenEll] Remark 1.5.1 —— **同型の spreading out**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★何を組み立てたか

`BaseChangeRatTower.lean` までで**射の側**は出ていた:

    X_ℚ ⟶ X′  は  X ×_ℤ ℤ[1/n!] ⟶ X′  を経由する

★原文が実際に使うのは**同型の側**である:

    X_ℚ ≅ X′_ℚ  ⟹  ∃n,  X ×_ℤ ℤ[1/n!] ≅ X′ ×_ℤ ℤ[1/n!]

★★組み立ては 4 段:

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | `e.hom` を降ろして `α_n : X_n ⟶ X′_n` | `exists_factor_baseChangeRatTower` ＋ `liftDescent` |
| 2 | `e.inv` を降ろして `β_n : X′_n ⟶ X_n` | 同上（左右を入れ替える） |
| 3 | `α ≫ β` と `𝟙` が生成ファイバーで一致 | ★`Spec ℤ[1/n!]` 成分は**引き戻しの普遍性だけ**で出る |
| 4 | ゆえに有限段で一致 | `descent_unique_baseChangeRatTower`（Stacks 01ZC 単射側） |

★★★**mathlib の欠落は無い**——段 4 が
`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType` そのものである。

## ★★★★★★★★配管——ここで 4 回落ちた（`tools/lean-idioms.md` に追記した）

`(baseChangeRatTowerDiagram f).obj n` と `pullback (overRatTowerDiagram.obj n).hom f` は
**defeq だが構文が違う**。混ぜると `pullback.hom_ext` が作る目標の中で
`Category.assoc` すら発火せず、

    Note: The target expression is not type-correct under the `instances` transparency level

が出る。★**対象・射・錐の脚を全部「生の引き戻しの形」に揃える**のが直し方であり、
そのために `bcObj` / `bcPt` / `bcMap` / `bcLeg` を置いた。
★★片方だけ揃えても駄目である——`bcDescHom` の余域だけ生にすると、
今度は域（`(D f).map h'` の余域）が合わなくなる。**全部揃えて初めて通る。**

## ★★★★逸脱（明示）

原文は「`X_ℚ ≅ X′_ℚ`」としか書かない。★形式化では
**その同型が `Spec ℚ` 上であること**を仮定 `he` として明示的に受ける
（`e.hom ≫ pullback.fst = pullback.fst`）。
★★原文の文脈では `X_ℚ`・`X′_ℚ` は `ℚ`-スキームであり同型も `ℚ`-同型なので、
**追加の制限ではない**。受け取り側（`Remark 1.5.1` の本体）に影響しない。

## ★残っている段

`Remark 1.5.1` を閉じるには、あと 2 つ:
★因子 `D` の降下、★★`Σ` の外での conductor の一致。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★生の引き戻しの形に揃える

★これらは `baseChangeRatTowerDiagram` / `baseChangeRatTowerCone` と **defeq** である。
別名を置くのは中身を変えるためではなく、**構文を 1 つに揃えて `rw` を通すため**である。 -/

/-- ★段 `n` の底変換 `X ×_ℤ ℤ[1/n!]`（生の引き戻しの形）。 -/
noncomputable abbrev bcObj {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℕᵒᵖ) :
    Scheme.{0} := pullback (overRatTowerDiagram.obj n).hom f

/-- ★極限の頂点 `X_ℚ`（生の引き戻しの形）。 -/
noncomputable abbrev bcPt {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) :
    Scheme.{0} := pullback (overRatTowerCone.pt).hom f

/-- ★段を下げる射。 -/
noncomputable def bcMap {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    {m n : ℕᵒᵖ} (h : m ⟶ n) : bcObj f m ⟶ bcObj f n :=
  (baseChangeRatTowerDiagram f).map h

/-- ★錐の脚 `X_ℚ ⟶ X_n`。 -/
noncomputable def bcLeg {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℕᵒᵖ) :
    bcPt f ⟶ bcObj f n := (baseChangeRatTowerCone f).π.app n

@[simp] theorem bcMap_snd {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    {m n : ℕᵒᵖ} (h : m ⟶ n) :
    bcMap f h ≫ pullback.snd (overRatTowerDiagram.obj n).hom f
      = pullback.snd (overRatTowerDiagram.obj m).hom f := baseChangeMap_snd f h

@[simp] theorem bcMap_fst {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    {m n : ℕᵒᵖ} (h : m ⟶ n) :
    bcMap f h ≫ pullback.fst (overRatTowerDiagram.obj n).hom f
      = pullback.fst (overRatTowerDiagram.obj m).hom f ≫ (overRatTowerDiagram.map h).left :=
  baseChangeMap_fst f h

theorem bcMap_comp {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    {l m n : ℕᵒᵖ} (h : l ⟶ m) (h' : m ⟶ n) :
    bcMap f (h ≫ h') = bcMap f h ≫ bcMap f h' := (baseChangeRatTowerDiagram f).map_comp h h'

@[simp] theorem bcLeg_snd {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℕᵒᵖ) :
    bcLeg f n ≫ pullback.snd (overRatTowerDiagram.obj n).hom f
      = pullback.snd (overRatTowerCone.pt).hom f := by
  show (baseChangeRatTowerCone f).π.app n ≫ _ = _
  simp only [baseChangeRatTowerCone, Functor.mapCone_π_app, Over.forget_map,
    Over.pullback_map_left]
  exact pullback.lift_snd _ _ _

@[simp] theorem bcLeg_fst {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℕᵒᵖ) :
    bcLeg f n ≫ pullback.fst (overRatTowerDiagram.obj n).hom f
      = pullback.fst (overRatTowerCone.pt).hom f ≫ (overRatTowerCone.π.app n).left := by
  show (baseChangeRatTowerCone f).π.app n ≫ _ = _
  simp only [baseChangeRatTowerCone, Functor.mapCone_π_app, Over.forget_map,
    Over.pullback_map_left]
  exact pullback.lift_fst _ _ _

/-- ★錐の脚は段を下げる射と両立する。 -/
theorem bcLeg_map {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    {m n : ℕᵒᵖ} (h : m ⟶ n) : bcLeg f m ≫ bcMap f h = bcLeg f n :=
  (baseChangeRatTowerCone f).w h

/-! ## ★★★★★★段を任意に下げられる降下射 -/

/-- ★★★★**降下した射を段 `n` で見たもの** ——
`X_i ⟶ X′` を段 `n ≤ i` へ引き戻し、`X′_n` へ持ち上げる。

★★これで「どの段でも同じ材料から作った」形になり、
段を下げたときの両立（`bcDescHom_map`）が言えるようになる。 -/
noncomputable def bcDescHom {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} (a : bcObj f i ⟶ X') {n : ℕᵒᵖ} (h : n ⟶ i) : bcObj f n ⟶ bcObj f' n :=
  liftDescent f f' (bcMap f h ≫ a)

@[simp] theorem bcDescHom_snd {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} (a : bcObj f i ⟶ X') {n : ℕᵒᵖ} (h : n ⟶ i) :
    bcDescHom f f' a h ≫ pullback.snd (overRatTowerDiagram.obj n).hom f'
      = bcMap f h ≫ a := liftDescent_snd _ _ _

@[simp] theorem bcDescHom_fst {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} (a : bcObj f i ⟶ X') {n : ℕᵒᵖ} (h : n ⟶ i) :
    bcDescHom f f' a h ≫ pullback.fst (overRatTowerDiagram.obj n).hom f'
      = pullback.fst (overRatTowerDiagram.obj n).hom f := liftDescent_fst _ _ _

/-- ★★★★★**段を下げても降下射は同じもの** ——
`X_m ⟶ X_n ⟶ X′_n` と `X_m ⟶ X′_m ⟶ X′_n` が一致する。

★これが無いと「段 `n` で示した等式を段 `m` へ運ぶ」ができない。 -/
theorem bcDescHom_map {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} (a : bcObj f i ⟶ X') {n m : ℕᵒᵖ} (h : n ⟶ i) (h' : m ⟶ n) :
    bcMap f h' ≫ bcDescHom f f' a h = bcDescHom f f' a (h' ≫ h) ≫ bcMap f' h' := by
  apply pullback.hom_ext
  · rw [Category.assoc, bcDescHom_fst, bcMap_fst, Category.assoc, bcMap_fst,
      ← Category.assoc, bcDescHom_fst]
  · rw [Category.assoc, bcDescHom_snd, Category.assoc, bcMap_snd, bcDescHom_snd,
      bcMap_comp, Category.assoc]

/-! ## ★★★★★★★★生成ファイバーでの一致 -/

/-- ★★★★★★**降下射は生成ファイバーで元の射に戻る**。

★仮定 `hfst` は「`u` が `Spec ℚ` 上の射である」こと、
`hsnd` は「`a` が `u` の降下である」こと。 -/
theorem bcLeg_bcDescHom {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} (a : bcObj f i ⟶ X')
    (u : bcPt f ⟶ bcPt f')
    (hfst : u ≫ pullback.fst (overRatTowerCone.pt).hom f'
      = pullback.fst (overRatTowerCone.pt).hom f)
    (hsnd : bcLeg f i ≫ a = u ≫ pullback.snd (overRatTowerCone.pt).hom f')
    {n : ℕᵒᵖ} (h : n ⟶ i) :
    bcLeg f n ≫ bcDescHom f f' a h = u ≫ bcLeg f' n := by
  apply pullback.hom_ext
  · rw [Category.assoc, bcDescHom_fst, bcLeg_fst, Category.assoc, bcLeg_fst,
      reassoc_of% hfst]
  · rw [Category.assoc, bcDescHom_snd, Category.assoc, bcLeg_snd,
      ← Category.assoc, bcLeg_map, hsnd]

/-- ★★★★★★★**降下射の合成が恒等になる判定** ——
段 `n` で「`X` 成分が一致する」ことが段 `m` へ引き戻して言えれば、
段 `m` では**合成そのものが恒等**である。

★機構は引き戻しの普遍性: `Spec ℤ[1/m!]` 成分は `bcDescHom_fst` から自動、
`X` 成分だけが仮定 `H` である。 -/
theorem bcDescHom_comp_id_of {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    {i j : ℕᵒᵖ} (a : bcObj f i ⟶ X') (b : bcObj f' j ⟶ X)
    {n m : ℕᵒᵖ} (h : n ⟶ i) (h2 : n ⟶ j) (h' : m ⟶ n)
    (H : bcMap f h' ≫ bcDescHom f f' a h ≫ bcDescHom f' f b h2 ≫
          pullback.snd (overRatTowerDiagram.obj n).hom f
        = bcMap f h' ≫ pullback.snd (overRatTowerDiagram.obj n).hom f) :
    bcDescHom f f' a (h' ≫ h) ≫ bcDescHom f' f b (h' ≫ h2) = 𝟙 (bcObj f m) := by
  apply pullback.hom_ext
  · rw [Category.assoc, bcDescHom_fst, bcDescHom_fst, Category.id_comp]
  · rw [Category.assoc, bcDescHom_snd, bcMap_comp, Category.assoc,
      ← bcDescHom_snd f' f b h2, ← Category.assoc, ← bcDescHom_map,
      Category.assoc, H, bcMap_snd, Category.id_comp]

/-! ## ★★★★★★★★★★同型の spreading out -/

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 同型の spreading out（両向きの射の形）**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`X_ℚ ≅ X′_ℚ`（`Spec ℚ` 上）なら、**有限段 `n` で `X ×_ℤ ℤ[1/n!] ≅ X′ ×_ℤ ℤ[1/n!]`**。

★★★組み立てはファイル冒頭の 4 段。★**mathlib の欠落は無い**。

★★仮定 `he` は「同型が `Spec ℚ` 上であること」——
原文の文脈（`X_ℚ`・`X′_ℚ` は `ℚ`-スキーム）では追加の制限ではない。

★★★★★添字を 1 つに揃える所で `ℕᵒᵖ` が**前順序**であること（hom が subsingleton）を使う
——`min` で 2 回落とした後、2 本の道 `n ⟶ i` が自動的に等しくなる。 -/
theorem exists_iso_baseChangeRatTower {X X' : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    [CompactSpace X'] [QuasiSeparatedSpace X']
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFinitePresentation f] [LocallyOfFinitePresentation f']
    (e : bcPt f ≅ bcPt f')
    (he : e.hom ≫ pullback.fst (overRatTowerCone.pt).hom f'
      = pullback.fst (overRatTowerCone.pt).hom f) :
    ∃ (n : ℕᵒᵖ) (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj f' n ⟶ bcObj f n),
      φ ≫ ψ = 𝟙 (bcObj f n) ∧ ψ ≫ φ = 𝟙 (bcObj f' n) := by
  classical
  have he' : e.inv ≫ pullback.fst (overRatTowerCone.pt).hom f
      = pullback.fst (overRatTowerCone.pt).hom f' := by
    rw [← he, ← Category.assoc, e.inv_hom_id, Category.id_comp]
  -- ★段 1: `e.hom` を降ろす
  obtain ⟨i, a, hai, -⟩ := exists_factor_baseChangeRatTower f f'
    (e.hom ≫ pullback.snd (overRatTowerCone.pt).hom f')
  -- ★段 2: `e.inv` を降ろす
  obtain ⟨j, b, hbj, -⟩ := exists_factor_baseChangeRatTower f' f
    (e.inv ≫ pullback.snd (overRatTowerCone.pt).hom f)
  have hai' : bcLeg f i ≫ a = e.hom ≫ pullback.snd (overRatTowerCone.pt).hom f' := hai
  have hbj' : bcLeg f' j ≫ b = e.inv ≫ pullback.snd (overRatTowerCone.pt).hom f := hbj
  -- ★共通の段へ落とす
  have hki : IsCofiltered.min i j ⟶ i := IsCofiltered.minToLeft i j
  have hkj : IsCofiltered.min i j ⟶ j := IsCofiltered.minToRight i j
  set k : ℕᵒᵖ := IsCofiltered.min i j
  have hα : ∀ {n : ℕᵒᵖ} (h : n ⟶ i), bcLeg f n ≫ bcDescHom f f' a h = e.hom ≫ bcLeg f' n :=
    fun h => bcLeg_bcDescHom f f' a e.hom he hai' h
  have hβ : ∀ {n : ℕᵒᵖ} (h : n ⟶ j), bcLeg f' n ≫ bcDescHom f' f b h = e.inv ≫ bcLeg f n :=
    fun h => bcLeg_bcDescHom f' f b e.inv he' hbj' h
  -- ★段 3: 生成ファイバーでは合成が恒等
  have keyAB : bcLeg f k ≫ bcDescHom f f' a hki ≫ bcDescHom f' f b hkj = bcLeg f k := by
    rw [← Category.assoc, hα hki, Category.assoc, hβ hkj, ← Category.assoc,
      e.hom_inv_id, Category.id_comp]
  have keyBA : bcLeg f' k ≫ bcDescHom f' f b hkj ≫ bcDescHom f f' a hki = bcLeg f' k := by
    rw [← Category.assoc, hβ hkj, Category.assoc, hα hki, ← Category.assoc,
      e.inv_hom_id, Category.id_comp]
  -- ★段 4: ゆえに有限段で一致（両向き）
  obtain ⟨mA, hA1, hA2, hmA⟩ := descent_unique_baseChangeRatTower f f
    (bcDescHom f f' a hki ≫ bcDescHom f' f b hkj ≫
      pullback.snd (overRatTowerDiagram.obj k).hom f)
    (specZIsTerminal.hom_ext _ _)
    (pullback.snd (overRatTowerDiagram.obj k).hom f)
    (specZIsTerminal.hom_ext _ _)
    ((reassoc_of% keyAB) (pullback.snd (overRatTowerDiagram.obj k).hom f))
  obtain ⟨mB, hB1, hB2, hmB⟩ := descent_unique_baseChangeRatTower f' f'
    (bcDescHom f' f b hkj ≫ bcDescHom f f' a hki ≫
      pullback.snd (overRatTowerDiagram.obj k).hom f')
    (specZIsTerminal.hom_ext _ _)
    (pullback.snd (overRatTowerDiagram.obj k).hom f')
    (specZIsTerminal.hom_ext _ _)
    ((reassoc_of% keyBA) (pullback.snd (overRatTowerDiagram.obj k).hom f'))
  -- ★★`ℕᵒᵖ` は前順序なので 2 本の道は等しい
  have hA12 : hA1 = hA2 := Subsingleton.elim _ _
  have hB12 : hB1 = hB2 := Subsingleton.elim _ _
  subst hA12; subst hB12
  have hAB : bcDescHom f f' a (hA1 ≫ hki) ≫ bcDescHom f' f b (hA1 ≫ hkj)
      = 𝟙 (bcObj f mA) := bcDescHom_comp_id_of f f' a b hki hkj hA1 hmA
  have hBA : bcDescHom f' f b (hB1 ≫ hkj) ≫ bcDescHom f f' a (hB1 ≫ hki)
      = 𝟙 (bcObj f' mB) := bcDescHom_comp_id_of f' f b a hkj hki hB1 hmB
  -- ★★★両者を共通の段 `n` へ運ぶ
  have hn1 : IsCofiltered.min mA mB ⟶ mA := IsCofiltered.minToLeft mA mB
  have hn2 : IsCofiltered.min mA mB ⟶ mB := IsCofiltered.minToRight mA mB
  set n : ℕᵒᵖ := IsCofiltered.min mA mB
  have hABn : bcDescHom f f' a (hn1 ≫ hA1 ≫ hki) ≫ bcDescHom f' f b (hn1 ≫ hA1 ≫ hkj)
      = 𝟙 (bcObj f n) :=
    bcDescHom_comp_id_of f f' a b (hA1 ≫ hki) (hA1 ≫ hkj) hn1
      (by rw [(reassoc_of% hAB) (pullback.snd (overRatTowerDiagram.obj mA).hom f)])
  have hBAn : bcDescHom f' f b (hn2 ≫ hB1 ≫ hkj) ≫ bcDescHom f f' a (hn2 ≫ hB1 ≫ hki)
      = 𝟙 (bcObj f' n) :=
    bcDescHom_comp_id_of f' f b a (hB1 ≫ hkj) (hB1 ≫ hki) hn2
      (by rw [(reassoc_of% hBA) (pullback.snd (overRatTowerDiagram.obj mB).hom f')])
  have e1 : hn2 ≫ hB1 ≫ hki = hn1 ≫ hA1 ≫ hki := Subsingleton.elim _ _
  have e2 : hn2 ≫ hB1 ≫ hkj = hn1 ≫ hA1 ≫ hkj := Subsingleton.elim _ _
  rw [e1, e2] at hBAn
  exact ⟨n, _, _, hABn, hBAn⟩

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 同型の spreading out（同型の形）**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

    X_ℚ ≅ X′_ℚ  ⟹  ∃n,  X ×_ℤ ℤ[1/n!] ≅ X′ ×_ℤ ℤ[1/n!]

★★これが原文が `Remark 1.5.1` の証明で使う spreading out の**同型の側**である。 -/
theorem exists_iso_of_generic_iso {X X' : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    [CompactSpace X'] [QuasiSeparatedSpace X']
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFinitePresentation f] [LocallyOfFinitePresentation f']
    (e : bcPt f ≅ bcPt f')
    (he : e.hom ≫ pullback.fst (overRatTowerCone.pt).hom f'
      = pullback.fst (overRatTowerCone.pt).hom f) :
    ∃ n : ℕᵒᵖ, Nonempty (bcObj f n ≅ bcObj f' n) := by
  obtain ⟨n, φ, ψ, h1, h2⟩ := exists_iso_baseChangeRatTower f f' e he
  exact ⟨n, ⟨⟨φ, ψ, h1, h2⟩⟩⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` には
★因子 `D` の降下、★★`Σ` の外での conductor の一致、が残っている。 -/

def exists_iso_baseChangeRatTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(同型の spreading out——因子 D の降下と conductor の一致は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_iso_of_generic_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(X_ℚ ≅ X′_ℚ ⟹ ∃n, X ×_ℤ ℤ[1/n!] ≅ X′ ×_ℤ ℤ[1/n!])",
    sectionId := "genell-rem-1-5-1" }

def exists_iso_of_generic_iso.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_factor_baseChangeRatTower(射の側の spreading out)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_factor_baseChangeRatTower") 9,
    .citation "[ABC3]" "descent_unique_baseChangeRatTower(降下の一意性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.descent_unique_baseChangeRatTower") 9,
    .citation "[mathlib]" "Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType(Stacks 01ZC 単射側)"
      (.inMathlib "AlgebraicGeometry.Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType") 9,
    .implicitStep
      ("★★原文は「X_ℚ ≅ X′_ℚ」としか書かない。形式化では同型が Spec ℚ 上であることを" ++
       "仮定 he として明示的に受けた。原文の文脈では X_ℚ・X′_ℚ は ℚ-スキームなので" ++
       "追加の制限ではない") 9,
    .implicitStep
      ("★★★★★添字を 1 つに揃える所で ℕᵒᵖ が前順序であること(hom が subsingleton)を使う。" ++
       "min で 2 回落とした後、2 本の道 n ⟶ i が自動的に等しくなる") 9,
    .implicitStep
      ("★残る段 1: 因子 D の降下——X の上の算術因子を X ×_ℤ ℤ[1/n!] へ降ろす") 9,
    .implicitStep
      ("★★残る段 2: Σ の外での conductor の一致。" ++
       "★これが揃えば Skeleton/GenEll/Section1.lean の remark_1_5_1 が受けている" ++
       "仮定 hagree を証明で置き換えられる") 9 ]

end ABC3.Found.GenEll
