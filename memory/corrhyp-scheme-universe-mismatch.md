---
name: corrhyp-scheme-universe-mismatch
description: HyperbolicCurveData.Space:Type(Type 0固定)にSchemeベースのデータを入れようとして「Application type mismatch」で衝突した。Schemeは宿命的にType 1に住む。
metadata:
  type: feedback
---

`Interface/CorrHyp/HyperbolicCurve.lean` の `HyperbolicCurveData.Space : Type`
(何も書かなければ `Type 0`)に、`AlgebraicGeometry.Scheme` ベースの具体的な
データ(`Over BaseK`、`BaseK : Scheme`)を入れようとして
`Application type mismatch ... of sort 'Type 1' but is expected to have type
?m.1 of sort 'Type'` で止まった(2026-09-04、[[corrhyp-track-goal]])。

**Why:** mathlib の `Scheme`(や `CommRingCat` 等、"大きな圏" の対象型)は
台となる `Type u` 値のデータ(層・環など)を**フィールドとして**含む構造体
なので、`Scheme.{u}` 自身は `Type (u+1)` に住む——これは避けられない
("大きな圏" の宿命であり、バグでも設計ミスでもない)。ゆえに
`Over (X : Scheme.{u}) : Type (u+1)` も同様。一方、interface 側の
`Space : Type` は暗黙に `Type 0` に固定されており、この2つは**原理的に
一致しえない**。

厄介なのは**エラーの出方**: 構造体リテラル `{ Space := Over BaseK, ... }`
の中で `Space := Over BaseK` 単独をチェックしてもエラーが出ないことがある
(unify が後続フィールドまで postpone されるらしく、実際に矛盾が表面化する
のは**別の**、`Space` を参照する後続フィールド(例: `ModuliStack := ...`)
の行になる)。★**「エラーが出た行」を鵜呑みにせず、`Space :=` の行を
単独で `#check (Over BaseK : Type)` のように孤立させて再現するのが早い**
——このとき初めて本当の場所(`Space` 自体)でエラーが出る。

**How to apply:**
- `Space`(や同様の「台」フィールド)に `Scheme`/`CommRingCat` 等の
  "大きな圏" のデータを載せる設計をするなら、**最初から** interface 側を
  universe 多相にしておく(`structure Foo.{u} where Space : Type u`)。
  既存の `FuchsianGroup`(`Type 0`)ベースの具体化には影響しない
  (`u := 0` で自動的にそのまま通る、`ULift` は要らない)。
- 逆に、後から気づいて直す場合も**恐れず直接 `Type` → `Type u` に書き換えて
  よい**——`structure` 自体を `.{u}` にするだけで済み、既存の具体化・
  既存の `variable (D : Foo)` を使う定理群には影響しない(Lean が
  universe を自動汎化する)。[[corrhyp-track-goal]] で実測済み
  (`lake build` で CorrHyp トラック全体を再ビルドし無傷を確認)。
- `Fuchsian`/`Aut` のように **大きくなくてよい**フィールドは、`Space` を
  上げても**そのまま `Type`(`Type 0`)に残してよい**——構造体の中で
  フィールドごとに universe を分けられる(全部を一律で上げる必要はない)。
