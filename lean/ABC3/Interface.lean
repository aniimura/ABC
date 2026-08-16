import ABC3.Meta.Claim
import ABC3.Interface.PGC.LocalFieldData
import ABC3.Interface.IUTchIII.PilotObjects
import ABC3.Interface.GenEll.ArithLineBundle
/-!
# Interface — まだ無い基礎の型

Track A(骨格)が上から必要とし、Track B(基礎)がまだ供給していないものを、
**`structure` として**置く場所。

## 規則

- **`axiom` を使わない**。`axiom` は検出不能な破滅、`structure` は検出可能な空虚
  (`ABC3/Meta/Calibration.lean` に実演がある)。
- ここに置く各 `structure X` は、次のどちらかを伴わなければならない
  (`tools/check.mjs` が G2 として検査する):
  - `X.nonvacuous : Nonempty X` — 仮説を満たす実例を実際に構成した
  - `X.waiting : ABC3.Meta.WaitingFor` — まだ構成できない。**何を待っているかを書く**
- 各フィールドは原典の逐語に裏付けられていなければならない(G1)。
  裏付けの無いフィールドを足すと、後で「自明化は原典のせいか我々のせいか」を
  切り分けられなくなる。
-/
