import ABC3.Interface.Arakelov.LineBundle

/-!
# 退化封じの検査 —— **`pullback` の公理だけでは層の引き戻しに決まらない**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★2026-08-18 に見つかった 2 つ目の穴

`Interface/Arakelov/LineBundle.lean` の `PicardData` には
`pullback` / `pullback_mul` / `pullback_id` / `pullback_comp` があり、
別に `sheafOf` / `sheafOf_one` / `sheafOf_mul` / `sheafOf_injective` /
`sheafOf_surjective` があった。

★★★**しかしどの条件も `pullback` と `sheafOf` を結んでいなかった。**

## ★★本ファイルが機械的に示すこと

引き戻しの 3 公理(`mul` / `id` / `comp`)だけを取り出した構造 `WeakPullback` は、

    Pic X := Multiplicative ℤ(すべての `X` で同じ)
    pullback f := id

という**幾何と一切関係のない witness** を持つ。
★したがって 3 公理は `pullback` を「層の引き戻し」に**決めない**。

★★★これが `PicardData.sheafOf_pullback` を追加した理由である。

## ★注意(限界の明示)

★本ファイルは「3 公理が `pullback` を決めない」ことを示すのであって、
「`PicardData` 全体が退化 witness を持つ」ことを示すのではない。
★★後者を機械的に示すには `PicardData` の witness が 1 つ要るが、
それはまだ無い(B1 が未完だから)。**この 1 点は検査ではなく設計判断である。**
-/

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory

/-- ★`PicardData` から引き戻しの 3 公理だけを取り出した構造。 -/
structure WeakPullback where
  /-- `Pic` の台。 -/
  Pic : Scheme.{0} → Type
  /-- 群構造。 -/
  group : (X : Scheme.{0}) → CommGroup (Pic X)
  /-- 引き戻し。 -/
  pullback : {X Y : Scheme.{0}} → (X ⟶ Y) → Pic Y → Pic X
  /-- 引き戻しは積を保つ。 -/
  pullback_mul : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : Pic Y),
    pullback f (@HMul.hMul _ _ _
        (@instHMul _ (group Y).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _ (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (pullback f L) (pullback f M)
  /-- 恒等射の引き戻しは恒等。 -/
  pullback_id : ∀ {X : Scheme.{0}} (L : Pic X), pullback (𝟙 X) L = L
  /-- 合成の引き戻しは引き戻しの合成。 -/
  pullback_comp : ∀ {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (L : Pic Z),
    pullback (f ≫ g) L = pullback f (pullback g L)

/-- ★★★★**幾何と無関係な witness** —— 定数関手 + 恒等引き戻し。

★★★これが 3 公理をすべて満たす以上、3 公理は `pullback` を
「層の引き戻し」に**決めていない**。 -/
def constWitness : WeakPullback where
  Pic _ := Multiplicative ℤ
  group _ := inferInstance
  pullback _ := id
  pullback_mul _ _ _ := rfl
  pullback_id _ := rfl
  pullback_comp _ _ _ := rfl

/-- ★★★★★**穴は実在した** —— 引き戻しの 3 公理は幾何と無関係な witness を通す。

★これが `PicardData.sheafOf_pullback`(2026-08-18 追加)の存在理由である。 -/
theorem weakPullback_nondegenerate_hole :
    ∃ W : WeakPullback, ∀ (X : Scheme.{0}) (f : X ⟶ X) (L : W.Pic X),
      W.pullback f L = L :=
  ⟨constWitness, fun _ _ _ => rfl⟩

/-! ## ★出典の紐付け(`.src`) -/

def weakPullback_nondegenerate_hole.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しの公理が層の引き戻しを決めないこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Check.Arakelov
