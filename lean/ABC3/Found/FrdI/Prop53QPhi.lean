/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Rlf
import ABC3.Found.FrdI.Prop53Base

/-!
# [FrdI] `ℚ·Φ^birat` を `Φ^pf` の中で立てる

原文 (FrdI p.103):
> if C is of Frobenius-isotropic type, then there is a natural 1-commutative

★★`Proposition 5.3` の第 2 文は `(𝒞^un-tr)^pf` の有理関数の単系を

  `ℚ·Φ^birat = Φ^birat ⊗_ℤ ℚ = (Φ^birat)^pf`

と **3 通りに**書く。★我々はすでに 1 つめを `sPhiBiratOn ℚ≥0`(`Prop53Rlf.lean`、
`ℚ≥0 ⊗_ℕ Φ` の中の ℚ-スパン)として持っているが、
**`𝒞^pf` の因子の単系は `Pf Φ`(商模型)である**(`pfRootPre`)。
そこで同じものを **`Pf Φ` の中で**立て直すのが本ファイルである。

## ★立て方

`Pf M` は**一意可除**なので、`ℚ`-スパンは**飽和**として書ける:

  `qPhiBiratOn G d := {x ∈ Gp(Pf Φ(d)) | ∃ k ≥ 1, k·x ∈ Φ^birat(d) の像}`

★これが部分群になるのは分母の**共通分母**を取るだけである
(`k·x ∈ S`, `l·y ∈ S` なら `(k·l)·(x+y) = l·(k·x) + k·(l·y) ∈ S`)。

★★★2 つの模型が一致すること(`qPhiBiratOn` と `sPhiBiratOn ℚ≥0` の対応)は
`pfEquivPfT`(`Def24PfT.lean`、`Pf M ≃+ ℚ≥0 ⊗_ℕ M`)で移すだけである。

## ★残り(記録)

原文の 3 つめ「`(Φ^birat)^pf`」との一致、すなわち
**`(Φ^pf)^birat = ℚ·Φ^birat`** はまだである。
★`ℚ` の出どころ(`𝒞^pf` の零因子の分母)は `Prop55PfDiv.lean` で測った。
残るのは双有理射の `Div^gp`(`sliceDivGpOf`)へその計算を移す段である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

section QPhi

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

variable (P) in
/-- ★`Φ^birat(d)` の `Pf` への像。 -/
noncomputable def phiBiratPfImage (G : Frobenioid P) (d : D) :
    AddSubgroup (Gp (Pf (Φ.val d))) :=
  AddSubgroup.map (gpMap _ (Pf.of (M := Φ.val d))) (phiBiratOn G d)

theorem mem_phiBiratPfImage {G : Frobenioid P} {d : D} {y : Gp (Φ.val d)}
    (hy : y ∈ phiBiratOn G d) :
    gpMap _ (Pf.of (M := Φ.val d)) y ∈ phiBiratPfImage P G d :=
  ⟨y, hy, rfl⟩

variable (P) in
/-- ★★★★**`ℚ·Φ^birat ⊆ (Φ^pf)^gp`** —— `Pf` の中で立てた版。

★`Pf M` は一意可除なので、`ℚ`-スパンは**「正の整数倍で `Φ^birat` の像に入る元の全体」**
(飽和)として書ける。 -/
noncomputable def qPhiBiratOn (G : Frobenioid P) (d : D) :
    AddSubgroup (Gp (Pf (Φ.val d))) where
  carrier := {x | ∃ k : ℕ+, ((k : ℕ+) : ℕ) • x ∈ phiBiratPfImage P G d}
  zero_mem' := ⟨1, by
    show ((1 : ℕ+) : ℕ) • (0 : Gp (Pf (Φ.val d))) ∈ _
    rw [smul_zero]
    exact (phiBiratPfImage P G d).zero_mem⟩
  add_mem' := by
    rintro x y ⟨k, hk⟩ ⟨l, hl⟩
    refine ⟨k * l, ?_⟩
    have h : ((k * l : ℕ+) : ℕ) • (x + y)
        = ((l : ℕ+) : ℕ) • (((k : ℕ+) : ℕ) • x) + ((k : ℕ+) : ℕ) • (((l : ℕ+) : ℕ) • y) := by
      rw [smul_add, smul_smul, smul_smul]
      push_cast
      rw [mul_comm ((l : ℕ+) : ℕ) ((k : ℕ+) : ℕ)]
    rw [h]
    exact (phiBiratPfImage P G d).add_mem
      ((phiBiratPfImage P G d).nsmul_mem hk _) ((phiBiratPfImage P G d).nsmul_mem hl _)
  neg_mem' := by
    rintro x ⟨k, hk⟩
    refine ⟨k, ?_⟩
    rw [smul_neg]
    exact (phiBiratPfImage P G d).neg_mem hk

theorem mem_qPhiBiratOn_iff {G : Frobenioid P} {d : D} {x : Gp (Pf (Φ.val d))} :
    x ∈ qPhiBiratOn P G d ↔ ∃ k : ℕ+, ((k : ℕ+) : ℕ) • x ∈ phiBiratPfImage P G d :=
  Iff.rfl

/-- ★★`Φ^birat` の像はそのまま `ℚ·Φ^birat` に入る(`k = 1`)。 -/
theorem phiBiratPfImage_le_qPhiBiratOn (G : Frobenioid P) (d : D) :
    phiBiratPfImage P G d ≤ qPhiBiratOn P G d := fun _ hx => ⟨1, by simpa using hx⟩

/-- ★★★**`ℚ·Φ^birat` は可除**(正の整数で割れる)——`ℚ`-スパンであることの中身。 -/
theorem qPhiBiratOn_of_nsmul_mem {G : Frobenioid P} {d : D} {x : Gp (Pf (Φ.val d))} (k : ℕ+)
    (h : ((k : ℕ+) : ℕ) • x ∈ qPhiBiratOn P G d) : x ∈ qPhiBiratOn P G d := by
  obtain ⟨l, hl⟩ := h
  refine ⟨l * k, ?_⟩
  have h2 : ((l * k : ℕ+) : ℕ) • x = ((l : ℕ+) : ℕ) • (((k : ℕ+) : ℕ) • x) := by
    rw [smul_smul]
    push_cast
    rfl
  rw [h2]
  exact hl

/-- ★★`𝒟` の射で引き戻る(部分関手性)。 -/
theorem qPhiBiratOn_map (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') {x : Gp (Pf (Φ.val d'))} (hx : x ∈ qPhiBiratOn P G d') :
    gpMap _ (Pf.map (Φ.map f)) x ∈ qPhiBiratOn P G d := by
  obtain ⟨k, y, hy, hxy⟩ := hx
  refine ⟨k, gpMap _ (Φ.map f) y, phiBiratOn_map G hiso hfn f hy, ?_⟩
  have h1 : gpMap _ (Pf.map (Φ.map f)) (((k : ℕ+) : ℕ) • x)
      = ((k : ℕ+) : ℕ) • gpMap _ (Pf.map (Φ.map f)) x := map_nsmul _ _ _
  rw [← h1, ← hxy, ← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply,
    ← gpMap_comp]
  refine congrArg (fun t : Φ.val d' →+ Pf (Φ.val d) => gpMap _ t y) ?_
  ext z
  show Pf.map (Φ.map f) (Pf.of z) = Pf.of (Φ.map f z)
  rfl

/-- ★`Pf.map` の逆(合成が恒等なら群化しても恒等)。 -/
theorem gpMap_pfMap_inv {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (f : M →+ N) (g : N →+ M) (h : g.comp f = AddMonoidHom.id M) (x : Gp (Pf M)) :
    gpMap _ (Pf.map g) (gpMap _ (Pf.map f) x) = x := by
  have hc : (Pf.map g).comp (Pf.map f) = AddMonoidHom.id (Pf M) := by
    ext z
    induction z using Pf.inductionOn with | _ m a =>
    show Pf.map g (Pf.map f (Pf.mk m a)) = Pf.mk m a
    rw [Pf.map_mk, Pf.map_mk]
    exact congrArg (fun t => Pf.mk t a) (congrArg (fun t : M →+ M => t m) h)
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, hc, gpMap_id]
  rfl

/-- ★★同型に沿っては両向き(`mem_sPhiBiratOn_iso` の `Pf` 版)。 -/
theorem mem_qPhiBiratOn_iso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') [IsIso f] (x : Gp (Pf (Φ.val d'))) :
    gpMap _ (Pf.map (Φ.map f)) x ∈ qPhiBiratOn P G d ↔ x ∈ qPhiBiratOn P G d' := by
  refine ⟨fun h => ?_, fun h => qPhiBiratOn_map G hiso hfn f h⟩
  have hcomp : (Φ.map (inv f)).comp (Φ.map f) = AddMonoidHom.id (Φ.val d') := by
    ext y
    show Φ.map (inv f) (Φ.map f y) = y
    rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id]
  have h2 := qPhiBiratOn_map G hiso hfn (inv f) h
  rwa [gpMap_pfMap_inv _ _ hcomp] at h2

end QPhi

/-! ## ★2. `sPhiBiratOn ℚ≥0`(テンソル模型)との一致 -/

section Bridge

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★単系の同型は群化の同型を与える(在庫に無かったので置く)。 -/
noncomputable def gpEquiv {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (e : M ≃+ N) :
    Gp M ≃+ Gp N where
  toFun := gpMap _ (e : M →+ N)
  invFun := gpMap _ (e.symm : N →+ M)
  left_inv x := by
    have hc : (e.symm : N →+ M).comp (e : M →+ N) = AddMonoidHom.id M := by
      ext y; exact e.symm_apply_apply y
    rw [← AddMonoidHom.comp_apply, ← gpMap_comp, hc, gpMap_id]
    rfl
  right_inv x := by
    have hc : (e : M →+ N).comp (e.symm : N →+ M) = AddMonoidHom.id N := by
      ext y; exact e.apply_symm_apply y
    rw [← AddMonoidHom.comp_apply, ← gpMap_comp, hc, gpMap_id]
    rfl
  map_add' x y := map_add _ x y

@[simp] theorem gpEquiv_apply {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (e : M ≃+ N)
    (x : Gp M) : gpEquiv e x = gpMap _ (e : M →+ N) x := rfl

/-- ★`Pf Φ(d) ≃+ ℚ≥0 ⊗_ℕ Φ(d)` を `Gp` へ持ち上げたもの。 -/
noncomputable def gpPfEquivPfT (M : Type w) [AddCommMonoid M] :
    Gp (Pf M) ≃+ Gp (PfT M) :=
  gpEquiv (pfEquivPfT (M := M))

end Bridge

/-! ## ★3. `ℚ·Φ^birat ⟶ ℝ·Φ^birat` —— `Proposition 5.3` の図式の右の縦の矢印の材料 -/

section ToRlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★自然数倍とスカラー倍は一致する(`Gp` の層で)。 -/
theorem sSmulGp_natCast {S : Type} [CommSemiring S] {N : Type w} [AddCommMonoid N]
    [Module S N] (n : ℕ) (z : Gp N) : sSmulGp ((n : S)) z = n • z := by
  obtain ⟨a, b, rfl⟩ : ∃ a b : N, z = toGp N a - toGp N b := by
    induction z using AddLocalization.induction_on with
    | H y => exact ⟨y.1, (y.2 : N), (eq_sub_of_add_eq (mk_add_toGp N y.1 y.2)).symm ▸ rfl⟩
  rw [sSmulGp_apply, map_sub, gpMap_toGp, gpMap_toGp]
  show toGp N ((n : S) • a) - toGp N ((n : S) • b) = n • (toGp N a - toGp N b)
  rw [Nat.cast_smul_eq_nsmul, Nat.cast_smul_eq_nsmul, smul_sub, toGp_nsmul, toGp_nsmul]

/-- ★★★**`S·Φ^birat` は可除**(`S` が半体なら) —— スカラー倍で閉じているから。 -/
theorem sPhiBiratOn_of_nsmul_mem {S : Type} [Semifield S] [CharZero S] (G : Frobenioid P) {d : D}
    {x : Gp (ScT S (Φ.val d))} (k : ℕ+)
    (h : ((k : ℕ+) : ℕ) • x ∈ sPhiBiratOn S G d) : x ∈ sPhiBiratOn S G d := by
  have hk : ((((k : ℕ+) : ℕ) : S)) ≠ 0 := by
    have : ((k : ℕ+) : ℕ) ≠ 0 := k.ne_zero
    exact_mod_cast Nat.cast_ne_zero.mpr this
  have h1 : sSmulGp ((((k : ℕ+) : ℕ) : S))⁻¹ (((k : ℕ+) : ℕ) • x) ∈ sPhiBiratOn S G d :=
    sPhiBiratOn_smul_mem G _ h
  rwa [← sSmulGp_natCast (S := S) ((k : ℕ+) : ℕ) x, sSmulGp_mul,
    inv_mul_cancel₀ hk, sSmulGp_one] at h1

/-- ★`Pf M ≃+ ℚ≥0 ⊗_ℕ M` は `Pf.of` を `toSc` へ移す。 -/
theorem pfEquivPfT_of {M : Type w} [AddCommMonoid M] (m : M) :
    pfEquivPfT (Pf.of m : Pf M) = toSc (S := NNRat) m := by
  show pnatInv (1 : ℕ+) ⊗ₜ[ℕ] m = (1 : NNRat) ⊗ₜ[ℕ] m
  rw [show pnatInv (1 : ℕ+) = (1 : NNRat) from by
    show (((1 : ℕ+) : ℕ) : NNRat)⁻¹ = 1
    norm_num]

/-- ★`Φ^pf = Pf Φ` から `Φ^rlf = ℝ≥0 ⊗_ℕ Φ` への標準の写像。 -/
noncomputable def pfToRlfHom (M : Type w) [AddCommMonoid M] : Pf M →+ ScT NNReal M :=
  (scBase (NNRat.castHom NNReal)).comp ((pfEquivPfT (M := M)) : Pf M →+ PfT M)

@[simp] theorem pfToRlfHom_of {M : Type w} [AddCommMonoid M] (m : M) :
    pfToRlfHom M (Pf.of m) = toSc m := by
  show scBase (NNRat.castHom NNReal) (pfEquivPfT (Pf.of m)) = _
  rw [pfEquivPfT_of, scBase_toSc]

/-- ★★★★★**`ℚ·Φ^birat` は `ℝ·Φ^birat` へ移る** ——
`Proposition 5.3` の図式の右の縦の矢印 `𝒞^pf ⟶ 𝒞^rlf` の単系の側。

★分母 `k` を払って `Φ^birat` の像に落とし、`ℝ≥0` 側では
スカラー倍で割り戻せる(`sPhiBiratOn_of_nsmul_mem`)。 -/
theorem qPhiBiratOn_map_rlf (G : Frobenioid P) {d : D} {x : Gp (Pf (Φ.val d))}
    (hx : x ∈ qPhiBiratOn P G d) :
    gpMap _ (pfToRlfHom (Φ.val d)) x ∈ sPhiBiratOn NNReal G d := by
  obtain ⟨k, y, hy, hyk⟩ := hx
  refine sPhiBiratOn_of_nsmul_mem G k ?_
  have h1 : ((k : ℕ+) : ℕ) • gpMap _ (pfToRlfHom (Φ.val d)) x
      = gpMap _ (pfToRlfHom (Φ.val d)) (((k : ℕ+) : ℕ) • x) := (map_nsmul _ _ _).symm
  rw [h1, ← hyk, ← AddMonoidHom.comp_apply, ← gpMap_comp]
  have hc : (pfToRlfHom (Φ.val d)).comp (Pf.of (M := Φ.val d))
      = (toSc : Φ.val d →+ ScT NNReal (Φ.val d)) := by
    ext z
    exact pfToRlfHom_of z
  rw [hc]
  exact mem_sPhiBiratOn_of_phiBiratOn G hy

end ToRlf

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.3` の `ℚ·Φ^birat` を `Φ^pf`(商模型)の中で立てたもの。 -/
def qPhiBiratOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — ℚ·Φ^birat(Φ^pf = Pf Φ の中の飽和として)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
