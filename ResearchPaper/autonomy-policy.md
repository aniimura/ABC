# 無人実行の線引き

2026-09-05 制定。`ResearchPaper/orchestration.md` の運用規則。

**「無人」= 人が画面を見ていない時間**を指す。人が居て判断を委ねられている場合は、
本体セッションがその場で決めてよい(2026-09-05 に実際そうした)。

## 0. 無人 `git push` は**許可**(2026-09-05、ユーザー判断)

理由:コミットはすべてゲート(全体ビルド + `check.mjs` + 文字化け)を通してから push しており、
`Found/` に `sorry` を残さない規約もあるので、壊れたものが master に載る経路が実質ない。
並行セッション(ABC3b/c)も同じ運用。

## 1. 無人でやってよいこと

- `node tools/frontier.mjs` で持ち場を選ぶ
- `lean-search` / `math-planner` を走らせる(どちらも書き込まない)
- **既に方針が決まっているノード**の実装(`pgc-goal.md` のノード表にあるもの)
- 対象モジュール単体のビルド
- 節目の全体ゲート(★**実装 agent が全員止まってから**。`check.mjs` は内部で
  `lake build` を走らせるのでロックで待たされる)
- ゲートが緑なら **commit + push**
- `meta-backlog.md` / `decisions-pending.md` / `memory/` への書き込み
- `tools/lean-idioms.md` への失敗形の追記

## 2. 必ず人を待つこと(★`decisions-pending.md` に積んで、止まらず次へ進む)

- **Skeleton の statement の修理**(退化の修理)——原典の読み替えを伴う
- **経路の選択**(A/B/C のような分岐)
- **メタ提案のマージ**(隔離 worktree からの取り込み)
- **`Interface` への場の追加・削除**——巻き添え範囲が広い
- **`CLAUDE.md` / `.mcp.json` / `.claude/settings.json` / `tools/hooks/` の変更**
- 原典の**逸脱**を新たに導入するとき

★**止まらないこと。** 判断が要るものに当たったら `decisions-pending.md` に積み、
`frontier.mjs` の次の startable へ進む。着手不可ノードを飛ばすのと同じ理屈。

## 3. 停止条件(これに当たったら波を止めて人を待つ)

- `node tools/check.mjs --brief` の NG が **13 件を超えた**
- `lake build ABC3` が**2 回連続で失敗**した(1 回目は並行セッションとの衝突を疑って再実行)
- 落ちたモジュールが**自分たちの下流**である
- `node tools/graph.mjs` のノード数が想定外に動いた

## 4. 上限

- 同時起動 **5 体**(`orchestration.md` §0)
- **main tree に書く agent は 1 波につき 1 体**まで(ビルドのロック競合を避ける)。
  残りは読み取り専用か隔離 worktree にする
- ★例外:互いに**別の新規ファイル**しか触らない実装 agent は同時に走らせてよい。
  ただし `Found.lean` / `Check.lean` への import 追加は**報告させて本体がまとめて行う**

## 5. 報告の義務

無人で回した波も、**何を判断し何を保留したか**を残す:

- commit メッセージに測定値(ジョブ数・NG 件数)を書く
- 保留した判断は `decisions-pending.md` に理由つきで
- メタの測定は `meta-backlog.md` に書き戻す
- 在庫の測定は `memory/` に(★「不在」の測定は**再実行できる形**で。
  2026-09-05 に 4 件が覆った)
