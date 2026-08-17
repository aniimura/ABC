---
name: lean-rebind-morphisms-clean-types
description: `F.map g` の型が `X ⟶ Y` に見えず `≫` も `IsIso` も落ちるとき、平の射へ束縛し直す。
metadata:
  type: feedback
---

**`(toBiratCat P G).map g` のような「関手の像」を `≫` で繋ぐと、型検査が
`instances` 透過度で落ちることがある。**

★原因: `(toBiratCat P G).obj A` や `(coaPreObj P G A).obj` は `A` と**定義的に等しい**が
**構文が違う**。`Quiver.Hom` の unification は `instances` 透過度で走るので、
そこを跨げない。

★★**巻き添えで `IsIso` のインスタンス探索まで落ちる** ——
`haveI : IsIso (F.map a)` を置いても直後の `inv (F.map a)` で
「failed to synthesize」になり、★**`IsIso (𝟙 X)` すら合成できなくなる**。
`maxHeartbeats` を上げても直らない(深さの問題ではない)。

## ★★★手当て(2026-08-17 に実測、[FrdI] `Prop 4.4` で 1 ターン溶かした)

1. ★**平の射へ束縛し直す**:
   `obtain ⟨aa, haa⟩ : ∃ aa : X ⟶ Y, aa = F.map g := ⟨_, rfl⟩`
   以後は `aa` だけを使い、必要なときだけ `rw [haa]` で戻す。
2. ★**`inv` を使わない** —— 逆射は `IsIso.out` で平の射として取り出し、
   打ち消しは `f ≫ (f' ≫ t) = t` の形の補題を手で用意する
   (`(Category.assoc _ _ _).symm.trans (congrArg (· ≫ t) h) |>.trans (Category.id_comp t)`)。
3. ★**`obtain ⟨Z, φ, rfl⟩` で射を代表元に潰さない** ——
   潰すと `HomBirat.mk Z φ` の型が `X ⟶ Y` に見えなくなり、`𝟙 X ≫ −` すら組めない。
   `hZφ : mk Z φ = f` のまま持ち、必要な所で `rw [← hZφ]`。

★同じ理由で `rw [Category.assoc]` / `rw [Category.comp_id]` が
「そこにあるのに当たらない」。★**`congrArg` と `calc` の項で書く**。

関連: [[lean-fullsubcat-procedures]]、[[lean-build-check-discipline]]。
