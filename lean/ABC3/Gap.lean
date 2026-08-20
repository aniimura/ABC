import ABC3.Meta.Claim
import ABC3.Gap.FrdI.Section1
import ABC3.Gap.FrdI.Section2
import ABC3.Gap.FrdI.Section5
import ABC3.Gap.FrdI.Section6
import ABC3.Gap.FrdI.Prop44
import ABC3.Gap.GenEll.BDDirection
/-!
# Gap — 飛躍(追加仮説として型に出す)

原典が確立していない段を、**主張を弱めずに**扱うための場所(G5)。

「証明できない」で止めず、不足を明示的な仮説として型に出す:

```lean
structure Gap_3_6 (d : Lemma35Data) : Prop where
  extra : (足りない性質を、原典の語彙で正確に書いたもの)

theorem prop_3_6 (d : Lemma35Data) (H : Gap_3_6 d) : Conclusion d := ...
```

各 `structure G` は `G.record : ABC3.Meta.GapRecord` を伴う
(`tools/check.mjs` が検査する)。`GapRecord.falsifier` は必須——
「何が起きればこの分類が覆るか」を書けないうちは、まだ ③ ではない。
-/
