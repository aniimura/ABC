/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ModelSquare

/-!
# [FrdI] Corollary 5.4 の縦の矢印の継ぎ目 —— `𝔽_Φ` の上で一致すれば model でも一致する

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★★`Cor54Compat.lean` の ★★ で書いた「要る 1 本」がこれである。

model Frobenioid の射は **4 成分 `(base, div, deg, u)`**、対象は **2 成分 `(base, cls)`**。
`ModelData.toElem`(`𝒞 ⥤ 𝔽_Φ`)が見るのは `base` と `div` と `deg` だけで、
★**落ちるのは対象の `cls` と射の `u` の 2 つ、すなわち有理関数の側だけ**である。

## ★★本ファイルが閉じること

| 定理 | 中身 |
|---|---|
| `ModelData.elemIsoBase` | ★`toElem` の後の自然同型から**底の同型**を取り出す |
| `ModelData.squareOfElemIso` | ★★★**`toElem` の後で一致 ＋ `(cls, u)` の帳尻 ⟹ 元の関手が同型** |

★`𝔽_Φ` の同型が `Div = 0`・`deg = 1` であることは在庫
(`ElemFrobCat.div_eq_zero_of_isIso` / `ElemFrobCat.isIso_iff`)を使う。
`IsDivisorial = IsPreDivisorial ∧ IsSharp` の **sharp** が可逆元を `0` に潰す。

## ★★★なぜ「等式」ではなく「同型」なのか(記録)

`Cor54Compat.lean` の旧い記述は継ぎ目を**関手の等式**

  `(unTr_modelFrobenioid Fc₁ G₁).functor ⋙ (untrModelHomOver …).functor`
    `= psiUnTr Ψ … ⋙ (unTr_modelFrobenioid Fc₂ G₂).functor`

として書いていたが、★原文が言う **1-commutative は「同型」**である。
`modelType_equiv` は `Classical.choice` で作った圏同値なので**等式は望むべくもなく**、
同型で十分かつ正しい。

★★`ModelSquare.lean` の `squareOfBaseU` は**底の同型 `e` を外から与える**形だったが、
`Corollary 5.4` の継ぎ目では `e` は `Corollary 4.11, (iv)` の 1-可換図式
(`cor_4_11_iv_square`)と `modelTypeEquiv_comp_toElem` の合成として
★**`𝔽_Φ` の上の自然同型の形でしか手に入らない**。本ファイルはその形を直接受ける。

★これで継ぎ目に残るのは **`(cls, u)` の帳尻 —— `uu` / `hc` / `hu` の 3 条**、
すなわち `Φ^birat` の輸送(`Cor54Birat.lean` の `phiBiratOn_transport_of_cor411`)
の側だけになる。`base` / `div` / `deg` の側は**もう要らない**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u3 v3

namespace ModelData

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}
  {E : Type u3} [Category.{v3} E]

/-! ## ★1. `toElem` の後の自然同型から底の同型を取り出す -/

/-- ★★**`toElem` の後の自然同型が与える底の同型**。

★`(K₁ ⋙ toElem).obj X` の底は `(K₁.obj X).base` そのもの(`toElem` は `⟨A.base⟩`)なので、
`α.hom.app X` の `base` 成分がそのまま `(K₁.obj X).base ⟶ (K₂.obj X).base` である。
★★逆射も `α.inv.app X` の `base` で取れるので、**`asIso` を経由せず**
構造体で直に組む(`.hom` が `rfl` になるようにするため)。 -/
def elemIsoBase {K₁ K₂ : E ⥤ Obj M} (α : K₁ ⋙ toElem ≅ K₂ ⋙ toElem)
    (X : E) : (K₁.obj X).base ≅ (K₂.obj X).base where
  hom := ElemFrobCat.Hom.base (α.hom.app X)
  inv := ElemFrobCat.Hom.base (α.inv.app X)
  hom_inv_id := congrArg ElemFrobCat.Hom.base (α.hom_inv_id_app X)
  inv_hom_id := congrArg ElemFrobCat.Hom.base (α.inv_hom_id_app X)

@[simp] theorem elemIsoBase_hom {K₁ K₂ : E ⥤ Obj M} (α : K₁ ⋙ toElem ≅ K₂ ⋙ toElem)
    (X : E) : (elemIsoBase α X).hom = ElemFrobCat.Hom.base (α.hom.app X) := rfl

/-! ## ★2. 継ぎ目の本体 -/

/-- ★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目**(道具の側)——

**`𝔽_Φ` へ落として一致し、落ちた 2 成分(`cls` と `u`)の帳尻が合えば、
model Frobenioid への 2 つの関手は同型である。**

★`base` / `div` / `deg` は仮定 `α` が全部与える(★同型なので `div = 0`・`deg = 1`)。
残る `cls` と `u` を `uu` / `hc` / `hu` で受け取る。 -/
noncomputable def squareOfElemIso (hdivOn : M.phi.IsDivisorialOn)
    (neg : ∀ d : D, M.bmon.val d → M.bmon.val d)
    (hneg : ∀ (d : D) (x : M.bmon.val d), neg d x + x = 0)
    {K₁ K₂ : E ⥤ Obj M} (α : K₁ ⋙ toElem ≅ K₂ ⋙ toElem)
    (uu : ∀ X : E, M.bmon.val (K₁.obj X).base)
    (hc : ∀ X : E, (K₁.obj X).cls
      = M.phi.gpMapOn (elemIsoBase α X).hom (K₂.obj X).cls + M.divB _ (uu X))
    (hu : ∀ {X Y : E} (f : X ⟶ Y),
      M.bmon.map (K₁.map f).base (uu Y) + (K₁.map f).u
        = M.bmon.map (elemIsoBase α X).hom (K₂.map f).u
          + (((K₂.map f).deg : ℕ+) : ℕ) • (uu X)) :
    K₁ ≅ K₂ :=
  NatIso.ofComponents (fun X => objIsoOfBaseU neg hneg (elemIsoBase α X) (uu X) (hc X))
    (fun {X Y} f => by
      haveI : ∀ Z : E, IsIso (α.hom.app Z) := fun Z =>
        inferInstanceAs (IsIso (α.app Z).hom)
      -- ★同型なので `𝔽_Φ` で落ちる 2 成分は `0` と `1`
      have hd0 : ∀ Z : E, ElemFrobCat.Hom.div (α.hom.app Z) = 0 := fun Z =>
        ElemFrobCat.div_eq_zero_of_isIso (fun A => (hdivOn A).2) _
      have hg1 : ∀ Z : E, ElemFrobCat.Hom.deg (α.hom.app Z) = 1 := fun Z =>
        ((ElemFrobCat.isIso_iff (α.hom.app Z)).mp inferInstance).2.2
      -- ★`α` の自然性を 3 成分に割る
      have hnat := α.hom.naturality f
      have hbase : (K₁.map f).base ≫ (elemIsoBase α Y).hom
          = (elemIsoBase α X).hom ≫ (K₂.map f).base :=
        congrArg ElemFrobCat.Hom.base hnat
      have hdeg : (K₁.map f).deg = (K₂.map f).deg := by
        have h : ElemFrobCat.Hom.deg (α.hom.app Y) * (K₁.map f).deg
            = (K₂.map f).deg * ElemFrobCat.Hom.deg (α.hom.app X) :=
          congrArg ElemFrobCat.Hom.deg hnat
        rw [hg1 X, hg1 Y, one_mul, mul_one] at h
        exact h
      have hdiv : (K₁.map f).div
          = M.phi.map (elemIsoBase α X).hom (K₂.map f).div := by
        have h := congrArg ElemFrobCat.Hom.div hnat
        simp only [ElemFrobCat.comp_div, hd0, hg1, map_zero, zero_add, PNat.one_coe,
          one_smul, smul_zero, add_zero] at h
        exact h
      refine Hom.ext ?_ ?_ ?_ ?_
      · exact hbase
      · show M.phi.map (K₁.map f).base (0 : M.phi.val _)
            + ((1 : ℕ+) : ℕ) • (K₁.map f).div
          = M.phi.map (elemIsoBase α X).hom (K₂.map f).div
            + (((K₂.map f).deg : ℕ+) : ℕ) • (0 : M.phi.val _)
        rw [map_zero, zero_add, smul_zero, add_zero, PNat.one_coe, one_smul]
        exact hdiv
      · show (1 : ℕ+) * (K₁.map f).deg = (K₂.map f).deg * (1 : ℕ+)
        rw [one_mul, mul_one]
        exact hdeg
      · show M.bmon.map (K₁.map f).base (uu Y) + ((1 : ℕ+) : ℕ) • (K₁.map f).u
          = M.bmon.map (elemIsoBase α X).hom (K₂.map f).u
            + (((K₂.map f).deg : ℕ+) : ℕ) • (uu X)
        rw [PNat.one_coe, one_smul]
        exact hu f)

end ModelData

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目
(★**条つき**: 落ちる 2 成分 `cls` / `u` の帳尻を仮定として受け取っている)。 -/
def ModelData.squareOfElemIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(𝔽_Φ の上の自然同型から model の自然同型へ)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
