---
name: mathlib-absent-by-wrong-name
description: 「mathlib に無い」の誤判定 6 件。原因は毎回「我々が付けたい名前で引いた」。下付き分岐群は Ideal.inertia として在った
metadata:
  type: reference
---

**mathlib の「不在」を、我々が付けたい名前で引いて判定してはいけない。**
2026-09-05〜06 に **6 件が覆った**。毎回同じ回路である。

| # | 「無い」と書いたもの | 実は在った名前 | 場所 |
|---|---|---|---|
| 1 | `ULift.field` | — | [[mathlib-frobenius-element-exists]] |
| 2 | Frobenius 元の存在 | — | 同上 |
| 3 | `Ẑ` | — | 同上 |
| 4 | `CompactSpace Gal` | — | 同上 |
| 5 | `IsArithFrobAt` | — | 同上 |
| 6 | ★**下付き分岐群 `G_i`（上付き↔下付き変換の前提）** | **`Ideal.inertia` / `AddSubgroup.inertia`** | 下記 |

## ★6 件目（2026-09-06、[pGC] §2 の価格測定）

`Skeleton/PGC/Section2.lean` の `prop_2_2.needs` が
`re:lowerNumbering|LowerNumbering|upperNumbering|UpperNumbering|herbrand|Herbrand → 0`
を根拠に「下付き番号付けも 0 件」と記録していた。**名前は無いが対象は在る**:

- `Ideal.inertia (G) (I) : Subgroup G` — `RingTheory/Ideal/Defs.lean:152`（`abbrev`）
- `AddSubgroup.inertia (I) (G) : Subgroup G` — `Algebra/Group/Subgroup/Basic.lean:1066`
- `AddSubgroup.mem_inertia : σ ∈ I.inertia G ↔ ∀ x, σ • x - x ∈ I` — **`@[simp]`**
- ★`AddSubgroup.subgroupOf_inertia : (I.inertia G).subgroupOf H = I.inertia H` — **`rfl`**
  ⇒ 原文（Serre / Milne）の「下付きは部分群への移行と両立する」が**そのまま在る**
- `AddSubgroup.inertia_map_subtype`、`Ideal.inertiaEquiv`（群同型に沿った移送）

すなわち `G_i(L/K) := Ideal.inertia (L ≃ₐ[K] L) (𝔪_L^(i+1))` で**定義は 1 行**。
★**本当に不在なのは上付き番号付けと Herbrand の変換だけ**である
（`RingTheory/Valuation/RamificationGroup.lean` が 4 宣言しかないのは 3 回測って正しい。
誤りは「だから下付きも無い」の側）。

**Why:** 我々は原典の用語（`lowerNumbering`）で索引を引くが、mathlib は
**別の数学的定式化**（作用が I を法として自明になる元の集合 = inertia）で入れていることがある。
名前の一致は在庫の判定基準にならない。

**How to apply:**

- ★**「無い」と書く前に、名前ではなく `statement` で引く。**
  `.cache/mathlib-index.txt` は statement 欄を持つので、**結論の形**（`Subgroup G`、
  `σ • x - x ∈ I` のような**定義の中身**）で grep する。
  ★名前欄と statement 欄の**両方**で引くこと（[[mathlib-index-nonascii-truncation]]）。
- ★**数学的に同値な別定式化を 2〜3 個考えてから引く。**
  分岐群なら「lowerNumbering」「ramificationGroup」に加えて
  「inertia」「stabilizer」「fixedPoints」「acts trivially mod I」。
- `node tools/absent-recheck.mjs --try '<正規表現>'` は**再実行できる形**で残すためのもの。
  `.needs` の `.absent` にはこの正規表現を必ず書く（規約）。★ただし
  **その正規表現が我々の語彙でしかないなら、0 件は「不在」の証拠にならない**。
- ★`.absent` の記録を消さずに**追記で訂正する**。何をどう測って外したかが次の較正になる。
