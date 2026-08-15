import ABC3.Found.FrdI.Prop19

/-!
# [FrdI] Proposition 1.6 —— Categorical Fiber Products

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.27–p.28(**400 dpi 目視確認 2026-08-15、実装のための再読**)。
§0 の `CFP` は物理 p.17(同じく目視)。

原文 (FrdI p.27):
> Proposition 1.6.

原文 (FrdI p.27):
> Let D′ be a connected, totally epimorphic category; D′ →D a functor that

## ★規模の測定(目視)

| 条 | 内容 | 主張の数 |
|---|---|---|
| (i) | `𝔽_Φ' ≅ 𝔽_Φ ×_𝒟 𝒟'` | 1 |
| (ii) | `𝒞' = 𝒞 ×_𝒟 𝒟'` が Frobenioid | 1(＋`Definition 1.3` の 21 条) |
| (iii) | 5 クラスが射影で判定できる | 5 |
| (iv) | base-iso について 3 クラス ＋ `𝒪^▷` の全単射 | 4 |
| (v) | 11 個の対象タイプ | 11 |
| (vi) | 3 タイプ(★**片向きのみ**) | 3 |

★**証明は 1 段落**(11 行)。21 条は「routine verification」の一言で片づけられている。

## ★★仕分け —— `Istr` とは**保証の出どころが違う**

`Proposition 1.9, (v)` では `Definition 1.3, (vii), (b)`(isotropic からの射の終域は
isotropic)が**閉性**を与え、前向きの構成が自動で `𝒞^istr` に留まった。
**CFP に閉性は無い。** 代わりに効くのは:

★**新しく作る対象・射は必ず base-isomorphism に沿って現れ、
`𝒞'` の base-isomorphism とは「`𝒟'` 成分が同型」という意味**なので、
`𝒟'` 成分を `𝟙` に取れる(新規対象)か、入力から与えられる(`plBkEquiv` など)。

★したがって **`G : 𝒟' ⥤ 𝒟` の充満性は要らない**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']

/-! ## ★`Φ` の `𝒟'` への制限

原文 (FrdI p.27):
> monoid obtained by restricting Φ to D′. Then:

★**`FSM ↦ FSM` の仮定はここでだけ使う** —— `Definition 1.1, (ii), (b)` を
`Φ'` について確かめるのに要る。`(ii), (a)` は `Φ` のものがそのまま降りる。 -/

/-- `Φ` を関手 `G : 𝒟' ⥤ 𝒟` に沿って制限した monoid。 -/
def MonoidOn.restrict (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    MonoidOn.{v3, u3, w} D' where
  functor := G.op ⋙ Φ.functor
  charInj α := Φ.charInj (G.map α)
  fsmIso α hα := Φ.fsmIso (G.map α) (hG α hα)

@[simp] theorem MonoidOn.restrict_val (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) (A : D') :
    (Φ.restrict G hG).val A = Φ.val (G.obj A) := rfl

theorem MonoidOn.restrict_map (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
    {A B : D'} (α : B ⟶ A) (x : Φ.val (G.obj A)) :
    (Φ.restrict G hG).map α x = Φ.map (G.map α) x := rfl

/-! ## ★`𝒞' = 𝒞 ×_𝒟 𝒟'` -/

variable {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

/-- **`𝒞' = 𝒞 ×_𝒟 𝒟'`**。 -/
abbrev CfpCat (P : PreFrobenioid C Φ) (G : D' ⥤ D) : Type _ := CFP P.proj G

/-- ★`𝒞'` の射の `𝒞` 成分。 -/
abbrev CfpCat.fst {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y) :
    X.obj.left ⟶ Y.obj.left := f.hom.left

/-- ★`𝒞'` の射の `𝒟'` 成分。 -/
abbrev CfpCat.snd {P : PreFrobenioid C Φ} {G : D' ⥤ D} {X Y : CfpCat P G} (f : X ⟶ Y) :
    X.obj.right ⟶ Y.obj.right := f.hom.right

/-- ★**`𝒞'` は totally epimorphic**。

原文 (FrdI p.28):
> the fact that D′ is a totally epimorphic category implies immediately that C′ is as

★原文は `𝒟'` の側しか挙げないが、**`𝒞` の側も要る**(射は対なので、
両成分がそれぞれ epi でなければならない)。`𝒞` は Frobenioid なので既に totally epimorphic。 -/
theorem cfp_totEpi (P : PreFrobenioid C Φ) (G : D' ⥤ D) (hD' : IsTotallyEpimorphic D') :
    IsTotallyEpimorphic (CfpCat P G) := by
  intro X Y f
  constructor
  intro Z g h hgh
  haveI h1 : Epi (CfpCat.fst f) := P.totEpiC _ _ _
  haveI h2 : Epi (CfpCat.snd f) := hD' _ _ _
  have e1 : CfpCat.fst f ≫ CfpCat.fst g = CfpCat.fst f ≫ CfpCat.fst h :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hgh
  have e2 : CfpCat.snd f ≫ CfpCat.snd g = CfpCat.snd f ≫ CfpCat.snd h :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hgh
  exact InducedCategory.hom_ext
    (CommaMorphism.ext ((cancel_epi (CfpCat.fst f)).mp e1) ((cancel_epi (CfpCat.snd f)).mp e2))

end ABC3.Found.FrdI
