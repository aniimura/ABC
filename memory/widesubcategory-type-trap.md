---
name: widesubcategory-type-trap
description: WideSubcategory / toElem を通ると対象が {obj := A}.obj の形になり、A.base と構文上は別物になる。射の構成は素の型を取る補題に出す。
metadata:
  type: feedback
---

`WideSubcategory (coaPreProp P)` や `P.toElem.obj A` を通ると、同じ対象が
`{ obj := A }.obj.base` / `((modelPre h).toElem.obj A).base` の形で現れる。
`A.base` と**定義上は等しいが構文上は別物**なので、

- `rw` がパターンを見つけられない（`Did not find an occurrence of the pattern`）
- `HAdd (Gp (Φ.val A.base)) (Gp (Φ.val ((P.toElem.obj A)).base))` の合成に失敗する
- `rw` の motive が型検査に落ちる

**How to apply:**

1. ★**射の構成は「素の型を取る補題」に出す**。圏同値の証明では、その補題を
   適用するだけにする（`exists_hom_under` / `exists_hom_over` / `exists_hom_div`
   が `Thm52Frob.lean` の実例）。補題の中は全部きれいな型なので `rw` が通る。
2. `leOfHom` などで取り出した値は
   `obtain ⟨x, rfl⟩ : ∃ x : 素の型, x = x₀ := ⟨x₀, rfl⟩` で**型を取り直す**。
   仮定も `have hx : きれいな形 := hx₀` で言い直す（`Eq` は defeq で通る）。
3. `IsIso.eq_inv_comp` のような instance 引数つきの補題は、`haveI` で
   instance を用意したうえで `@` で明示的に渡す（`IsIso` は `Prop` なので
   証明無関係、別の instance を渡しても定義的に一致する）。

**Why:** 2026-08-18、`Theorem 5.2, (ii)` の (iii)(d) 圏同値でこれに 4 回連続で
つまずいた。原因は毎回同じで、**症状（`rw` 失敗・instance 合成失敗・motive 不整合）が
違って見えるだけ**だった。

関連: [[lean-build-check-discipline]]
