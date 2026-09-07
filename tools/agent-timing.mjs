#!/usr/bin/env node
// agent-timing.mjs —— M6(セッションの実挙動が観測できない)の観測点。
//
// ★発見(メタ第 23 回): agent の完了通知 `<usage><subagent_tokens>/<tool_uses>/<duration_ms></usage>` は
//   会話ログ `~/.claude/projects/<slug>/<session>.jsonl` に **そのまま残っている**。
//   さらに `<session>/subagents/agent-*.jsonl` に **子 agent の全 tool 呼び出しと brief** が残っている。
//   ⇒ 本体が通知から手で書き写す必要は無い。この道具が読み出す。
//
// 使い方:
//   node tools/agent-timing.mjs --list                  … 通知を新しい順に出す
//   node tools/agent-timing.mjs --stats                 … 1 ノードあたりの時間の分布
//   node tools/agent-timing.mjs --explain               … 何が時間を説明するか(族は固定)
//   node tools/agent-timing.mjs --estimate              … brief の見積 行数 と 実測 行数
//   node tools/agent-timing.mjs --cost                  … 自己申告 COST[…] を実測に突き合わせる
//   node tools/agent-timing.mjs --main                  … ★本体セッションの費用(道具ごと / 失敗形ごと / ファイル塊)
//   node tools/agent-timing.mjs --denominator           … ★★数えているのは「実行」か「書かれた字面」か(M142)
//   node tools/agent-timing.mjs --m149                  … ★★★事前登録(0b)の実行。★件数が届くまで検定しない
//   node tools/agent-timing.mjs --solo                  … ★★M159 「他と混ざらない命令だけ」を**規則を焼いて**数える
//   node tools/agent-timing.mjs --concurrency           … ★★★M166 同時実行数の上限 2 は妥当か(事前登録つき)
//   node tools/agent-timing.mjs --mcp-watch             … ★★M173 「MCP の名指し」の前後を見張る(★分母つき)
//   node tools/agent-timing.mjs --supply 2026-09-07     … ★M179 その日の持ち場が frontier に載っていたか
//   node tools/agent-timing.mjs --supply 2026-09-07 --history … ★★歴史の木を立てて測り直す(重い)
//   node tools/agent-timing.mjs --selftest              … 自己検査
//   共通: --type lean-prover  --day 2026-09-07  --name <正規表現>  --limit N  --json
//         --since 2026-09-07T16:35:29Z  --until …  … ★時刻で絞る(ISO。辞書順 = 時刻順)
//         --ms 919271,1086819,…  … duration_ms を並べて過去の表をそのまま再現する
//
// ★統計の作法(M90 に倣う):
//   - 比較の族 FAMILY を **データを見る前に**コードへ固定する。
//   - 判定は「言える / 言えない」だけ。★率や相関の**向きを断定しない**。
//   - 多重比較(Holm)で補正する。
//   - 交絡を必ず出す。標本が独立でないとき(束)は ICC / DEFF を出す。

// ★★★自己申告の書式(メタ第 24 回。M96 が「今の記録では測れない」と書いた穴を塞ぐ)
//
//     COST[<持ち場>]: <安|並|高>  — <一言>
//     COST[<鍵>]: <安|並|高> | 持ち場=<agent に配った呼び名>  — <一言>     ★鍵が違うときの橋
//
//   - `<持ち場>` は **agent に配ったときの呼び名**(`Y19e` `B2+B3` `Λ12` など)。
//     ★`GUESS/VERDICT`(autonomy-policy §4.7)と**同じ字面**にすること。ここだけが対応の鍵。
//   - ★★**枝番を付けない**(`LT-a` ではなく `LT`)。費用は**持ち場 1 つに 1 行**である
//     (見当は 1 持ち場に何本あってもよいが、かかった時間は 1 つしかない)。
//   - ★★**実測は「鍵が持ち場名と一致すること」に依存する。**★第 24 回の実測:
//     いま `decisions-pending.md` にある鍵 13 個のうち、**agent の呼び名に当たるのは 2 個だけ**
//     (`段1` `C49` `L9` `TL` `L12` `HS` `AE` `GC` `LT` `IR` `SR` は当たらない。略号だから)。
//     ⇒ 鍵を揃えられないときは `| 持ち場=局所Tate双対性の道` と**明示の橋**を書くこと。
//     ⇒ ★当たらなかった申告は `--cost` が**名指しで**出す(黙って落とさない)。
//   - 意味は「**見積に対して**思ったより安かった / 見積どおり / 思ったより高かった」。
//     ★`安` は「行数が少なかった」ではない(それは `--estimate` が既に測っている)。
//   - **報告が返ったときに 1 行**書く(`VERDICT` と同じ間合い)。書かなくてもゲートは落ちない。
//   - 書き場所は `ResearchPaper/decisions-pending.md`(`GUESS/VERDICT` と同じ本)。
//
// ★なぜ語の検索ではだめか(M96 の実測): 「安い」系の語を含む節は 695 節中 21 節あるのに、
//   測れた 21 件の持ち場に**紐づくのは 1 件だけ**だった。★節は持ち場に対応していない。
//   ⇒ **鍵付きの 1 行**でなければ分母が立たない。
//
// ★なぜ `unverified.mjs` ではなく**ここ**に置いたか(3 つ)
//   1. ★M96 が測れなかったのは「申告と**実測**の一致」である。実測(所要時間 / 行数 / 見積)を
//      持っているのは**この道具だけ**。unverified 側に置くと「申告が N 件ある」までしか言えず、
//      穴は塞がらない。
//   2. ★unverified.mjs は **Bonferroni の族を固定してある**(M90)。同じ本の中に別種の欄を
//      足すと m の意味が濁り、**既に出ている判定が動く**。★族を守るために分ける。
//   3. ★unverified.mjs は import すると即座に集計が走る作り(main の番人が無い)。
//      読み手として使えないので、写すか作り替えるかになる。★どちらも今回の目的より高い。
//   ⇒ ★**書式は GUESS/VERDICT と揃え、置き場所も同じ本にし、数える道具だけ分ける。**

import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import readline from 'node:readline';
import { execFileSync } from 'node:child_process';

// ════════════════════════════════════════════════════════════════════
// 0. 比較の族 —— ★データを見る前に固定する(本体の brief が挙げた候補そのもの)
// ════════════════════════════════════════════════════════════════════
export const FAMILY = [
  // ★★族の正典は `lines`(いまの wc -l)。`linesWritten`(当時の Write 本文)は**族に入れない**。
  //   ── メタ第 25 回が測って決めた。理由は ρ の大小ではない(★下記のとおり ρ では決められない)。
  //   (1) ★ρ では決められない: n=191 で lines ρ=0.260(n=131) / linesWritten ρ=0.263(n=81)。
  //       ★メタ第 24 回の「0.758 対 0.511」(n=50)の差は **n を増やしたら消えた**。
  //       ⇒ ρ を根拠に選ぶのは「結果を見てから選ぶ」ことでもあり、そもそも差が無い。
  //   (2) ★どちらが何件間違うかで決めた。両方取れる 81 件のうち値が違うのは 53 件。
  //       ・そのうち **40 件は その agent 自身が最後の Write の後に Edit した**ぶん。
  //         ⇒ ここでは `linesWritten` のほうが**間違い**(自分の仕事を数え落とす)。
  //       ・自分は Write 以後さわっていないのにずれるのは 13 件。うち 9 件は **2 行以下**の差。
  //         ⇒ 「別の agent に後から書き換えられた」で実質的に効くのは **4 件 / 81 = 5%**。
  //       ⇒ 誤差 5% の `lines` を採り、誤差 49% の `linesWritten` は採らない。
  //   (3) ★取りこぼしも lines が有利: `linesWritten` は Write を使わなかった agent で欠測し、
  //       n が 131 → 81 に落ちる。★さらに族を差し替えると Holm の階段が早く止まり、
  //       **cores と checkFails の判定が「言える」→「言えない」に落ちる**(副作用が大きい)。
  //   (4) ★本当の正解はどちらでもない:「完了時点の行数」= 最後の Write 本文 + その後の自分の
  //       Edit を当てた行数。★これは識別もでき取りこぼしも無い。★未実装(meta-backlog M111)。
  //       ★実装したら**測ってから**差し替えること。今の判定を根拠に先に差し替えないこと。
  { key: 'lines',      label: '行数(wc -l)' },
  { key: 'toolUses',   label: 'tool_uses' },
  { key: 'tokens',     label: 'subagent_tokens' },
  { key: 'cores',      label: '抽象核の本数' },
  { key: 'leanChecks', label: 'lean_check の回数' },
  { key: 'checkFails', label: 'lean_check が失敗した回数(≒ 一発で通らなかった)' },
  // ★★事前登録(メタ第 23 回 M96 が「族の外」と自己申告し、★次の波で確かめると書いた欄)。
  //   ★第 24 回が**データを見る前に**ここへ移した。★向きは印字しない。
  { key: 'estMid',     label: '見積中点(brief の「見積 N–M 行」)' },
];

/** ★欄ごとの最小件数。これ未満なら検定しない(★データを見る前に決める)。 */
export const MIN_N = 8;

// ════════════════════════════════════════════════════════════════════
// 0b. ★★★M149 の事前登録(メタ第 31 回 2026-09-08 01:4x JST。★データを見る前に書いた)
// ════════════════════════════════════════════════════════════════════
/**
 * ★背景。メタ第 30 回 M145 が測った: `lines`/`cores` が見る「代表ファイル」は
 *   「`.lean` の path が tool_use の入力に現れた回数が最大の本」なので **Read しかしていない本が選ばれる**。
 *   実測 **138 件中 47 件(34%)** が「自分では 1 度も Write/Edit していない本」だった。
 *   分母を `written` に直すと `cores` の Holm が 0.0049 → 0.4071 に落ちる(= 判定が消える)。
 *   ★だが v1 は**データを見た後に**決めた分母なので、その p は事前登録の p と同じ資格を持たない
 *   (メタ第 24 回 M102 と同じ「2 度目の覗き」)。
 *
 * ★★以下は M149 が指示した手順そのものである。★**この節を書いた時点で新しい標本は 1 件も見ていない。**
 *
 * ── (A) 定義の決定(★データに依らない。いま決める) ────────────────────────────
 *   `lines`/`cores` の分母は「**その agent 自身が Write/Edit した `.lean` のうち最も多く触った本**」
 *   (= `covariatesFromText(..., {fileSource:'written'})`)が**正しい**。
 *   ★理由は相関の大小ではない: 「行数」と名付けた共変量が **その agent が 1 行も書いていない本の行数**を
 *   指しているのは**測定の誤り**であって、モデルの選択ではない。★だから検定の結果に依らず決まる。
 *   ★★ただし **この波では既定を差し替えない**(M149 の手順 1 は「コメントとして事前登録する」)。
 *   差し替えは (B) の結果を添えて本体が決める。★`FAMILY` は 7 本のまま、`m` は動かさない。
 *
 * ── (B) 検定の対象(★これだけが経験的な問い) ──────────────────────────────
 *   問い: 「**抽象核の本数は所要時間を説明する**」(M101/M110 の判定)は、
 *         分母を (A) に直しても立つか。
 *   標本: ★**M145 が見ていない agent だけ**。= `ts > M149_CUTOFF` の通知。
 *   統計: `familyLadder`(既存・事前登録済み)。★族は 7 本のまま。
 *         Spearman ρ / 並べ替え p(seed 固定)/ Holm。判定は `Holm < 0.05` の 2 値のみ。向きは断定しない。
 *   一次: `cores`(v1)。★他の 6 欄と v0 は**参考**で、判定の根拠にしない。
 *
 * ── (C) 締切 T の決め方(★結果に触れずに機械で決めた) ─────────────────────────
 *   `M149_CUTOFF` = 本体チェックアウトの `tools/agent-timing.mjs`(1,960 行 = M146 を本体が採用した版)と
 *   `ResearchPaper/meta-backlog.md`(M145–M151 を書き足した版)の **mtime の遅いほう**を秒に切り上げた時刻。
 *   実測(メタ第 31 回の起動時): 2026-09-08 01:35:28.365 / .458 JST = 2026-09-07T16:35:28.4Z。
 *   ⇒ ★この時刻は **M145 が測った瞬間より必ず後**である(測ってから報告し、本体が採用したのだから)。
 *     ⇒ これより後に始まった agent を M145 が見たことは**ありえない**。★安全側に倒してある。
 *
 * ── (D) 停止規則(★これがいちばん大事。★逐次の覗きを禁じる) ─────────────────────
 *   ★標本は時間とともに増える。★n=10 で覗いてから n=40 で検定するのは**逐次の覗き**で、
 *   第一種の過誤が膨らむ。⇒ ★**n が `M149_NEED` に届くまで階段を計算しない**。
 *   `--m149` は届いていなければ **件数と不足数だけ**を印字して止まる(★道具が規則を強制する)。
 *
 *   `M149_NEED` = 134。導出(★データを使っていない):
 *     Fisher-z 近似 n = 3 + ((z_{α/2} + z_β)/atanh ρ)^2、両側、検出力 0.80、
 *     α = 0.05/7(Holm の最悪の段 = Bonferroni)、
 *     ρ = **0.30**(★「実務で意味のある最小の効果」として**データを見ずに**決めた。
 *       これ未満の ρ なら「核の本数で工数を見積もる」という使い道が立たない)。
 *     ⇒ z=2.6901 / 0.8416 → n = 133.2 ⇒ **134**。
 *     ★参考(同じ式): ρ=0.20 → 307 / ρ=0.25 → 195 / ρ=0.35 → 97 / ρ=0.40 → 73。
 *   ★★危険側を先に書く: 標本は独立でない(M95 の束 / ICC)。DEFF > 1 なので
 *     **134 は下限**である。★束が効いていれば実効 n はこれより小さい。
 *
 * ── (E) M149 の手順 3(M111 も同時に決める) ─────────────────────────────
 *   ★「完了時点の行数」(M111)は **メタ第 26 回 M117 が既に測って決着している**:
 *     M111 の誤り 26%(21/82) > `lines` の汚染 5%(M110 の 4/81)、しかも n が 137 → 82 に落ちる。
 *   ⇒ ★**M111 は族に入れないし正典にもしない。この波でも次の波でも再開しない。**
 *     ★行数の候補 3 つ(`lines` / `linesWritten` / M111)のうち残る論点は
 *     **分母(touched か written か)だけ**であり、それが上の (A)(B) である。
 *   ⇒ ★これで「3 つの候補を別々の波で試して m を増やす」ことを避けられる(M149 の手順 3 の趣旨)。
 *
 * ── (F) 実行不能のときの答え方 ────────────────────────────────────────
 *   ★n が届かないときは「あと何件」を印字して**止まる**。★それが正しい結果である。
 */
export const M149_CUTOFF = '2026-09-07T16:35:29Z';   // ★(C) で決めた締切。★書き換えたら selftest が鳴る
export const M149_NEED = 134;                        // ★(D) で決めた必要件数。★書き換えたら selftest が鳴る
export const M149_PRIMARY = 'cores';                 // ★一次の欄
export const M149_DENOM = 'written';                 // ★(A) で決めた分母

/** ★M145 が見ていない通知だけを返す純関数(★ts は ISO 文字列。辞書順 = 時刻順)。 */
export function freshRows(recs, cutoff = M149_CUTOFF) {
  return recs.filter(r => typeof r.ts === 'string' && r.ts > cutoff);
}

/**
 * ★`--since` / `--until` の絞り込み(★純関数にした理由: CLI に埋めると selftest が届かず、
 *   ★第 31 回の突然変異 A10「--since を無視する」が**素通りした**)。
 * ★境界: `since` は**含まない**(freshRows と同じ)、`until` は**含む**。
 */
export function filterByTime(recs, since = null, until = null) {
  let out = recs;
  if (since) out = out.filter(r => typeof r.ts === 'string' && r.ts > since);
  if (until) out = out.filter(r => typeof r.ts === 'string' && r.ts <= until);
  return out;
}

/**
 * ★M149 の停止規則を**道具として**強制する純関数。
 * rows は既に `fileSource:'written'` で共変量を付けたもの。
 * 返り値: { n, need, open, ladder }  ★open が false なら ladder は **null**(計算すらしない)。
 */
export function m149Gate(rows, need = M149_NEED, key = M149_PRIMARY) {
  const n = usableRows(rows, key).length;
  const open = n >= need;
  return { n, need, open, short: Math.max(0, need - n), ladder: open ? familyLadder(rows) : null };
}

// ══════════════════════════════════════════════════════════════════════════
// ★★★M160 —— M149 の標本が**貯まっているか**を見張る(メタ第 32 回。持ち場 4)
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★★**いちばん大きな危険は保存期間である**(第 31 回の名指し)。
 *   M149 は「あと 134 件貯まるまで検定しない」という停止規則で待っている。
 *   ★ところが ★**古い側のログが刈られると、待っている間に標本が減りうる。**
 *   ★そうなると `--m149` は永久に開かないまま、しかも**それに気づけない**
 *   (`--m149` は「あと N 件」としか言わないので、N が増えても減っても同じ顔をする)。
 *
 * ⇒ ★**件数の履歴を外に持ち、単調に増えているかを見る。**
 *
 * ★★**履歴の置き場所が肝**: 改善係は毎回まっさらな隔離 worktree で起動し、
 *   `.cache/` は `meta-setup` の整列対象外(SKIP_SEG)なので**次の起動に残らない**。
 *   ⇒ ★`ResearchPaper/m149-watch.json` に置く(`COPY_ROOTS` に入っているので整列で運ばれる)。
 *
 * ★見る量は 3 つ。★どれが動いても意味が違う:
 *   - `usable` … 検定に使える件数。★**減ったら刈られている**(赤)。
 *   - `afterT` … 締切 T より後の通知。★減ったら刈られている(赤)。
 *   - `first`  … 記録全体のいちばん古い時刻。★**進んだら古い側が消えている**(黄)。
 *     ★これは `usable` が減る**前に**出る早期の合図である(T より前が先に刈られるため)。
 */
export const M149_WATCH_REL = 'ResearchPaper/m149-watch.json';

/** ★M171: これ未満の間隔では `perDay` / `etaDays` を出さない(0.25 日 = 6 時間)。
 *  ★書き換えたら selftest が鳴る。 */
export const M149_WATCH_MIN_DAYS = 0.25;

/** ★純関数。いまの観測 1 件を作る。 */
export function m149Observation(recs, usable, at = new Date().toISOString()) {
  const ts = recs.map((r) => r.ts).filter((x) => typeof x === 'string').sort();
  return {
    at,
    nRecs: recs.length,
    afterT: freshRows(recs).length,
    usable,
    first: ts[0] ?? null,
    last: ts[ts.length - 1] ?? null,
  };
}

/**
 * ★純関数。前回と今回を比べて判定を返す。
 * 返り値 { level: 'ok'|'warn'|'alarm', lines: string[], perDay: number|null, etaDays: number|null }
 */
export function m149WatchVerdict(prev, cur, need = M149_NEED) {
  const lines = [];
  let level = 'ok';
  if (!prev) {
    return { level: 'first', lines: ['★初回の観測。次に叩いたときから増減が見える。'], perDay: null, etaDays: null };
  }
  const bump = (l) => { if (l === 'alarm' || (l === 'warn' && level === 'ok')) level = l; };
  if (cur.usable < prev.usable) {
    bump('alarm');
    lines.push(`★★★使える件数が **減った**(${prev.usable} → ${cur.usable})。★ログが刈られている。`);
  }
  if (cur.afterT < prev.afterT) {
    bump('alarm');
    lines.push(`★★★T より後の件数が **減った**(${prev.afterT} → ${cur.afterT})。★ログが刈られている。`);
  }
  if (prev.first && cur.first && cur.first > prev.first) {
    bump('warn');
    lines.push(`★★いちばん古い記録が **進んだ**(${prev.first} → ${cur.first})。`
      + '★古い側から消えている ⇒ ★このまま待つと T より後にも届く。');
  }
  const dtDays = (Date.parse(cur.at) - Date.parse(prev.at)) / 86400000;
  let perDay = null, etaDays = null;
  /* ★★M171(メタ第 34 回) —— **最小間隔の守り**。
   *   ★M170 の実測: 2 つの観測の間隔が **14.2 分**しかないのに
   *   `usable` が 1 → 2 と 1 件増えただけで `perDay = 101 件/日 ⇒ あと 1.3 日` が出た。
   *   ★外挿の分母が小さすぎる。★しかも**楽観に振れる**ので「まだ余裕がある」と読ませる ——
   *     M160 の目的(刈られていないかの見張り)に対して**危険側**である。
   *   ⇒ `dtDays < M149_WATCH_MIN_DAYS` なら perDay / etaDays を **出さない**(null)。
   *   ★`level` の判定(alarm / warn)は間隔に依らないので**そのまま**である。 */
  let tooShort = null;
  if (dtDays > 0 && dtDays < M149_WATCH_MIN_DAYS) {
    /* ★★この行は **異常の報せではない**(見張りの結論は別に出す)。
     *   ⇒ `lines` に直接積むと「★単調に増えている。刈られた形跡は無い。」を押し出してしまう。
     *   ★実測(2026-09-08、実データ 54.7 分): 最初の実装がまさにその行を消した。 */
    tooShort = `★間隔が短すぎる(${(dtDays * 24 * 60).toFixed(1)} 分 < ${(M149_WATCH_MIN_DAYS * 24).toFixed(0)} 時間)。`
      + '★速さ(件/日)と残り日数は**出さない**。★半日以上あけてから叩くこと。';
  } else if (dtDays > 0) {
    perDay = (cur.usable - prev.usable) / dtDays;
    const short = Math.max(0, need - cur.usable);
    if (short === 0) etaDays = 0;
    else if (perDay > 0) etaDays = short / perDay;
    else { etaDays = Infinity; bump('warn'); lines.push('★増えていない(この間隔では 0 件/日)。★届く見込みが立たない。'); }
  }
  if (level === 'ok' && !lines.length) lines.push('★単調に増えている。刈られた形跡は無い。');
  if (tooShort) lines.push(tooShort);          // ★M171: 見張りの結論の**後ろ**に添える
  return { level, lines, perDay, etaDays };
}

// ══════════════════════════════════════════════════════════════════════════
// ★★★M166 —— 「同時実行数の上限 2 は妥当か」の事前登録(メタ第 33 回 2026-09-08 02:5x JST)
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★★★この節は **outcome の値を 1 つも見る前に**書いた。
 *   見たのは (i) 子 agent の jsonl に timestamp があること、
 *   (ii) tool の**名前**の一覧(Bash 10608 / lean_check 2089 / lean_start 102 / …)だけである。
 *   ★error の字面・所要時間・同時数のどれも見ていない。
 *
 * ── 問い ────────────────────────────────────────────────────────────
 *   D27 の「実装 agent は最大 2」は測って決めた数ではない。★3〜4 にしてよいか。
 *
 * ── ★測れないことを先に書く(交絡) ──────────────────────────────────
 *   (1) ★**難しい波ほど同時本数が多い**。⇒「同時数 → 悪化」の**因果は測れない**。
 *       測れるのは相関と、下の (a)(b)(c) の**件数**だけ。
 *   (2) ★**露出の作り方そのものが所要時間と機械的に絡む**: 長く走る agent ほど
 *       誰かと重なりやすいので `overlapMax` は duration と自動的に相関する。
 *       ⇒ ★**一次の露出は `overlapAtStart`**(自分が**始まった瞬間**に走っていた本数)にする。
 *         これは配る側が実際に決めている量であり、自分の duration では動かない。
 *         `overlapMax` / `overlapMean` は**参考**。
 *   (3) ★標本は独立でない(同じ波の agent は束)。ICC / DEFF を必ず出す。
 *
 * ── 露出(exposure) ────────────────────────────────────────────────
 *   agent i の区間は [s_i, e_i]、e_i = Date.parse(ts)、s_i = e_i − durationMs。
 *   `overlapAtStart(i)` … s_i の瞬間に走っていた本数(★i を含む。1 なら「1 人きり」)
 *   `overlapMax(i)`     … 区間中の同時本数の最大(i を含む)
 *   `overlapMean(i)`    … 区間中の同時本数の時間平均(i を含む)
 *   ★母集団は 2 つ作る:
 *     `all`   … 全 agent(読み取り専用も含む)
 *     `impl`  … ★`agentType === 'lean-prover'` だけ。★**D27 が上限を掛けているのはこちら**。
 *   ★一次は **`impl` の `overlapAtStart`**。
 *
 * ── 転帰(outcome)—— ★族は 4 本。★後から欄を足さない ─────────────────
 *   (a) `mcpInfra`   … MCP(`mcp__abc3-lean__*`)の tool_result のうち **配管の失敗**の件数。
 *                      ★判定は下の `MCP_INFRA_RE`(★データを見ずに書いた字面)。
 *                      ★Lean の型エラーは**入らない**(それは仕事であって取り合いではない)。
 *   (b) `sharedFile` … その agent が Write/Edit した本のうち、**区間が重なる別の agent も
 *                      Write/Edit した本**の数(distinct)。★衝突の**機会**であって衝突そのものではない。
 *   (c) `lakeWaitMs` … `LAKE_RE` に当たる Bash 呼び出しの**実時間の中央値**
 *                      (tool_use の timestamp → tool_result の timestamp)。
 *   (d) `msPerTool`  … `durationMs / toolUses`。★`durationMs` そのものは (2) で機械的に絡むので一次にしない。
 *
 * ── 検定 ──────────────────────────────────────────────────────────
 *   Spearman ρ / 並べ替え p(seed 固定、20000 回)/ Holm(m = 4)。
 *   ★欄ごとの最小件数は `MIN_N`(8)。届かない欄は**検定しない**。
 *   ★判定は「言える / 言えない」の 2 値だけ。★**向きを断定しない**(M90 の規律)。
 *   ★★どの欄が「言える」になっても、それは**相関**であって「3 本にしてよい」の根拠にはならない。
 *     根拠になるのは (a)(b)(c) の**件数が水準ごとにどうなっているか**の表である。
 *
 * ── ★判断の規則(★データを見る前に決める。これが結論の出し方) ──────────
 *   `impl` の `overlapAtStart` が 3 以上の区画に **agent が MIN_N 件以上**あり、かつ
 *   その区画の (a)+(b)+(c) の**発火件数が 0 件**なら「★3 でも配管は壊れていない」と**書いてよい**。
 *   1 件でも出たら件数を名指しし、★**「2 のままにすべき」と書く**。
 *   ★水準 3 以上の agent が MIN_N に届かないなら「★**言えない。あと何件**」と書いて止まる。
 */
export const FAMILY_C = [
  { key: 'mcpInfra',   label: '(a) MCP 配管の失敗の件数' },
  { key: 'sharedFile', label: '(b) 重なる agent と同じ本を書いた数' },
  { key: 'lakeWaitMs', label: '(c) lake build の実時間(中央値ms)' },
  { key: 'msPerTool',  label: '(d) 1 tool あたりの実時間(ms)' },
];

/**
 * ★★★v1 —— **事前登録した**(a) の検出器。★★実データでの感度は **0/19 だった**(下記)。
 * ★捨てずに残す。理由: 「事前登録した通りに走らせたら何が出たか」を後から再現できるようにするため。
 * ★★この定数を (a) の判定に使ってはいけない(`MCP_INFRA_RE` を使うこと)。
 */
export const MCP_INFRA_RE_V1 = new RegExp(
  [
    'no (running|active)( lean)? (session|instance|server)',
    '(session|instance|environment) (not found|expired|closed|died|invalid|mismatch)',
    'MCP error',
    'connection (closed|refused|reset|lost)',
    'ECONNRE', 'EPIPE', 'ETIMEDOUT',
    'server (not connected|disconnected|crashed|died|unavailable)',
    'not connected',
    'lean_start .{0,40}(fail|error)',
    'imports? (do not match|mismatch)',
    '環境が(違|異な)', 'インスタンスが(違|無|な)', 'セッションが(切|無|な)',
    'imports が(違|合わ)',
    '再起動が必要', 'restart(ing)? the (lean|server|session)',
  ].join('|'), 'i');

/**
 * ★★★v2(★**事後**。データを見て直した。★確証ではなく探索である)——
 *   メタ第 33 回が v1 を実データに当てたら **2,307 件中 0 件**しか鳴らなかった。
 *   ★ところが同じ木に `エラー: REPL は処理中(直列にしか使えない)` が **19 件**ある。
 *   ★これは #236 の取り合いそのものであり、v1 は **感度がゼロ**だった。
 *   ⇒ ★**「0 件」は「起きていない」ではなく「検出器が英語しか知らなかった」**である。
 *   ★★教訓: **鳴らないことを確かめていない検出器の 0 は、数字ではない。**
 *     ⇒ `MCP_INFRA_FIXTURES` を置き、selftest が**実物の字面で鳴ることを毎回確かめる**。
 *
 * ★正常形は入れない: `OK (0.07 秒)` / `エラー N 件` (= Lean の型エラー。仕事であって取り合いではない)
 *   / `REPL を落とした(再生用の控えも捨てた)` (= `lean_reset` の成功) は**数えない**。
 */
export const MCP_INFRA_RE = new RegExp(
  [
    // ★実データから(メタ第 33 回に全数を数えて拾った 3 形)
    'REPL は処理中',                       // ★取り合いそのもの(直列にしか使えない)
    '秒で応答が無いので REPL を落とした',   // ★時間切れで落ちた(240〜600 秒)
    'まだ lean_start を呼んでいない',       // ★環境が消えている / 手順違い
    // ★v1 の英語側(将来の実装で出うる形。★いまの木では 1 件も鳴らない)
    'no (running|active)( lean)? (session|instance|server)',
    '(session|instance|environment) (not found|expired|closed|died|invalid|mismatch)',
    'MCP error',
    'connection (closed|refused|reset|lost)',
    'ECONNRE', 'EPIPE', 'ETIMEDOUT',
    'server (not connected|disconnected|crashed|died|unavailable)',
  ].join('|'), 'i');

/**
 * ★★検出器の**感度**を毎回確かめる実物の字面(メタ第 33 回に木から全数で拾った)。
 * ★`want: true` は鳴らなければならない、`false` は鳴ってはいけない。
 * ★★この表があるかぎり「(a) が 0 件」は**検出器が生きている上での 0** である。
 */
export const MCP_INFRA_FIXTURES = [
  { want: true,  n: 19, s: 'エラー: REPL は処理中(直列にしか使えない)' },
  { want: true,  n: 17, s: 'エラー: 600 秒で応答が無いので REPL を落とした。次の呼び出しで自動的に再起動して import を読み直す。' },
  { want: true,  n: 6,  s: '[{"type":"text","text":"まだ lean_start を呼んでいない。imports を指定して lean_start を呼ぶこと。"}]' },
  { want: false, n: 4,  s: '[{"type":"text","text":"REPL を落とした(再生用の控えも捨てた)。"}]' },
  { want: false, n: 0,  s: '[{"type":"text","text":"OK (0.07 秒)"}]' },
  { want: false, n: 0,  s: '[{"type":"text","text":"エラー 1 件 (1.01 秒)\\n\\nerror 15:53\\nunsolved goals"}]' },
  { want: false, n: 1,  s: 'Error: result (75,857 characters across 971 lines) exceeds maximum allowed tokens.' },
];

/** ★`lake build` の待ちとみなす Bash(★データを見る前に固定)。 */
export const LAKE_RE = /\blake\s+(build|env)\b|tools[\\/]build\.mjs/;

/** ★書き込み衝突の機会を数える対象(★データを見る前に固定)。 */
export const CONFLICT_RE = /\.lean$|lean-idioms\.md$/i;

/** ★区画(水準)の切り方。★データを見る前に固定。4 以上は 1 区画にまとめる。 */
export const CONC_BINS = [1, 2, 3, 4];
export const concBin = (k) => (k >= 4 ? 4 : k);

/**
 * ★純関数。区間の配列から同時本数の 3 つの露出を作る。
 * iv: [{ id, s, e }]  ★s <= e。s === e(duration 0)も落とさない。
 * 返り値: Map(id -> { atStart, max, mean })
 */
export function overlapStats(iv) {
  const out = new Map();
  for (const a of iv) {
    // 開始の瞬間に走っている本数(★自分を含む。境界は「開いていれば走っている」= s_j <= s_a <= e_j)
    let atStart = 0;
    for (const b of iv) if (b.s <= a.s && a.s <= b.e) atStart++;
    // 区間中の同時本数の最大と時間平均。★変化点は他の区間の端点だけ。
    const pts = [a.s, a.e];
    for (const b of iv) {
      if (b.e < a.s || b.s > a.e) continue;
      if (b.s > a.s && b.s < a.e) pts.push(b.s);
      if (b.e > a.s && b.e < a.e) pts.push(b.e);
    }
    pts.sort((x, y) => x - y);
    let max = 0, acc = 0;
    const span = a.e - a.s;
    for (let i = 0; i + 1 < pts.length; i++) {
      const mid = (pts[i] + pts[i + 1]) / 2;
      let k = 0;
      for (const b of iv) if (b.s <= mid && mid <= b.e) k++;
      if (k > max) max = k;
      acc += k * (pts[i + 1] - pts[i]);
    }
    if (pts.length < 2) { max = atStart; acc = 0; }
    out.set(a.id, { atStart, max: Math.max(max, atStart), mean: span > 0 ? acc / span : atStart });
  }
  return out;
}

/**
 * ★純関数。子 agent の jsonl 本文から M166 の転帰を作る。
 * 返り値 { mcpInfra, mcpCalls, lakeCalls, lakeWaitMs, lakeMsList, wrote:Set, leanStarts }
 */
export function concurrencyOutcomesFromText(text) {
  const c = { mcpInfra: 0, mcpCalls: 0, lakeCalls: 0, lakeWaitMs: NaN, lakeMsList: [],
              wrote: new Set(), leanStarts: 0, mcpV1: 0, kind: { busy: 0, timeout: 0, nostart: 0, other: 0 } };
  const pendingMcp = new Set();
  const pendingLake = new Map();   // tool_use_id -> t0(ms)
  for (const ln of text.split('\n')) {
    if (!ln) continue;
    let o; try { o = JSON.parse(ln); } catch { continue; }
    const msg = o.message; if (!msg) continue;
    const t = o.timestamp ? Date.parse(o.timestamp) : NaN;
    const content = Array.isArray(msg.content) ? msg.content : [];
    for (const b of content) {
      if (b.type === 'tool_use') {
        const n = b.name || '';
        if (n.startsWith('mcp__abc3-lean__')) { c.mcpCalls++; pendingMcp.add(b.id); }
        if (n.endsWith('lean_start')) c.leanStarts++;
        if (n === 'Bash' && typeof b.input?.command === 'string' && LAKE_RE.test(b.input.command)) {
          c.lakeCalls++;
          if (Number.isFinite(t)) pendingLake.set(b.id, t);
        }
        const fp = b.input && (b.input.file_path || b.input.filePath);
        if (typeof fp === 'string' && CONFLICT_RE.test(fp)
            && (n === 'Write' || n === 'Edit' || n === 'MultiEdit')) {
          c.wrote.add(fp.replace(/\\/g, '/').toLowerCase());
        }
      } else if (b.type === 'tool_result') {
        if (pendingMcp.has(b.tool_use_id)) {
          pendingMcp.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          if (MCP_INFRA_RE_V1.test(s)) c.mcpV1++;          // ★事前登録した検出器(比較のためだけに残す)
          if (MCP_INFRA_RE.test(s)) {
            c.mcpInfra++;
            if (/REPL は処理中/.test(s)) c.kind.busy++;
            else if (/秒で応答が無いので REPL を落とした/.test(s)) c.kind.timeout++;
            else if (/まだ lean_start を呼んでいない/.test(s)) c.kind.nostart++;
            else c.kind.other++;
          }
        }
        if (pendingLake.has(b.tool_use_id)) {
          const t0 = pendingLake.get(b.tool_use_id);
          pendingLake.delete(b.tool_use_id);
          if (Number.isFinite(t) && t >= t0) c.lakeMsList.push(t - t0);
        }
      }
    }
  }
  if (c.lakeMsList.length) c.lakeWaitMs = quantile(c.lakeMsList, 0.5);
  return c;
}

/**
 * ★★事前登録(メタ第 25 回)—— 自己申告 COST を FAMILY に足してよいか。
 * ★n=5 の時点では「言えない」ですらなく **判定を出さない**(--cost が既にそう振る舞う)。
 * ★足す条件を**データを見る前に**ここへ固定する:
 *   (a) 突き合わせが MIN_N(8)件以上、かつ
 *   (b) ★水準が偏っていないこと —— 各水準に 3 件以上。
 * ★(b) が要る理由: 2026-09-07 時点の 5 件は 安:4 / 並:1 / 高:0 で、
 *   件数だけ 8 に達しても「全部 安」なら順位相関は定義できず、**分散ゼロの欄**になる。
 * ★足すときは m が 7 → 8 になり、既存の判定が Holm で罰される。★それを承知で 1 度だけ足すこと。
 */
export function costFamilyReady(levels) {
  const ls = Array.isArray(levels) ? levels : [];
  if (ls.length < MIN_N) return false;
  return COST_LEVELS.every(l => ls.filter(x => x === l).length >= 3);
}

/**
 * 族の 1 欄について「検定に使える行」を返す。
 * ★★欄ごとの complete-case(listwise ではない)。理由: `estMid` は brief に見積を
 *   書かなかった持ち場で欠測する。listwise だと **1 件の欠測で欄ごと落ちる**ので、
 *   事前登録した欄が「測れなかった」に化けてしまう。
 * ★欠測が欄によって違うので、**n を欄ごとに印字する**(そうしないと比較が読めない)。
 */
export function usableRows(rows, key) {
  return rows.filter((r) => Number.isFinite(r[key]));
}

// ════════════════════════════════════════════════════════════════════
// 1. 純関数 —— 通知の解析
// ════════════════════════════════════════════════════════════════════

/** task-notification の本文 1 件を解析する。usage が無ければ null。 */
export function parseNotification(content) {
  if (typeof content !== 'string') return null;
  const g = (re) => re.exec(content)?.[1];
  const d = g(/<duration_ms>(\d+)<\/duration_ms>/);
  if (d === undefined) return null;
  const summary = g(/<summary>([\s\S]*?)<\/summary>/) ?? '';
  return {
    taskId: g(/<task-id>([^<]*)<\/task-id>/) ?? '',
    toolUseId: g(/<tool-use-id>([^<]*)<\/tool-use-id>/) ?? '',
    status: g(/<status>([^<]*)<\/status>/) ?? '',
    summary,
    // `Agent "…" finished` の中身。二重引用符は content 側で素のまま入っている。
    name: /^Agent "([\s\S]*)" finished\s*$/.exec(summary.trim())?.[1] ?? summary.trim(),
    durationMs: Number(d),
    toolUses: Number(g(/<tool_uses>(\d+)<\/tool_uses>/) ?? NaN),
    tokens: Number(g(/<subagent_tokens>(\d+)<\/subagent_tokens>/) ?? NaN),
  };
}

/** brief(子 agent の最初の user message)から見積 行数 を読む。 */
export function parseEstimate(brief) {
  if (typeof brief !== 'string') return null;
  // 「見積 **400–700 行**」「見積 120–250 行（Y4）」「見積 **約 300 行**」等
  const m = /見積[^0-9]{0,12}(\d{2,4})\s*[-–—~〜]\s*(\d{2,4})\s*行/.exec(brief);
  if (m) return { lo: +m[1], hi: +m[2] };
  const s = /見積[^0-9]{0,12}(\d{2,4})\s*行/.exec(brief);
  if (s) return { lo: +s[1], hi: +s[1] };
  return null;
}

/** .lean の本文から「抽象核」の節を数える(節見出し `/-! ## … 抽象核 …`)。 */
export function countAbstractCores(src) {
  if (typeof src !== 'string') return 0;
  let n = 0;
  for (const line of src.split(/\r?\n/)) {
    if (/^\s*\/-!\s*#+\s.*抽象核/.test(line)) n++;
  }
  return n;
}

/** `COST[<持ち場>]: 安` を 1 行から読む。書式の**見本**(鍵に <> を含む)は数えない。 */
export const COST_RE = /^\s*(?:[★☆*\-\s]*)COST\[([^\]]+)\]\s*[:：]\s*(.*)$/;
export const COST_LEVELS = ['安', '並', '高'];
export function parseCostLine(line) {
  const m = COST_RE.exec(String(line ?? ''));
  if (!m) return null;
  const key = m[1].trim();
  if (/[<>]/.test(key)) return null;              // ★書式の見本を数えない(M71/unverified と同じ穴)
  const lv = /^\s*(安|並|高)(?![^\s—\-|]) ?/.exec(m[2]) ?? /^\s*(安|並|高)/.exec(m[2]);
  if (!lv) return null;
  const rest = m[2].slice(lv[0].length);
  // ★明示の橋: `| 持ち場=<agent の呼び名>`(GUESS の `確度=` `検算=` と同じ流儀)
  const slot = /持ち場\s*=\s*([^|—]+)/.exec(rest)?.[1]?.trim() || null;
  return {
    key,
    level: lv[1],
    slot,                                  // ★null なら key で突き合わせる
    note: rest.replace(/\|?\s*持ち場\s*=\s*[^|—]+/, '').replace(/^\s*[—\-|]\s*/, '').trim(),
  };
}

/** 記録(複数)から COST を集める。★同じ鍵は後の記述で上書きする。 */
export function readCosts(paths) {
  const out = new Map(); const missing = [];
  for (const p of paths) {
    if (!fs.existsSync(p)) { missing.push(p); continue; }
    const lines = fs.readFileSync(p, 'utf8').split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      const c = parseCostLine(lines[i]);
      if (c) out.set(c.key, { ...c, file: p, line: i + 1 });
    }
  }
  return { costs: out, missing };
}

/**
 * 申告の鍵と agent の持ち場名を突き合わせる。
 * ★前方一致だけだと `Y19` が `Y19e …` を拾ってしまうので、**鍵の直後は英数字でない**ことを要る。
 */
export function matchesKey(description, key) {
  const d = String(description ?? '').trim();
  const k = String(key ?? '').trim();
  if (!k) return false;
  if (d === k) return true;
  if (!d.startsWith(k)) return false;
  return !/[0-9A-Za-z]/.test(d.slice(k.length, k.length + 1));
}

// ════════════════════════════════════════════════════════════════════
// 2. 純関数 —— 統計
// ════════════════════════════════════════════════════════════════════

/** type-7 分位点(numpy / R の既定)。xs は昇順でなくてよい。 */
export function quantile(xs, q) {
  const a = [...xs].sort((x, y) => x - y);
  if (a.length === 0) return NaN;
  if (a.length === 1) return a[0];
  const h = (a.length - 1) * q;
  const lo = Math.floor(h), hi = Math.ceil(h);
  return a[lo] + (h - lo) * (a[hi] - a[lo]);
}

/** 平均順位(同順位は平均)。 */
export function ranks(xs) {
  const idx = xs.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
  const r = new Array(xs.length);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const avg = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) r[idx[k][1]] = avg;
    i = j + 1;
  }
  return r;
}

export function pearson(x, y) {
  const n = x.length;
  if (n < 2) return NaN;
  let mx = 0, my = 0;
  for (let i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
  mx /= n; my /= n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    const a = x[i] - mx, b = y[i] - my;
    sxy += a * b; sxx += a * a; syy += b * b;
  }
  if (sxx === 0 || syy === 0) return NaN;
  return sxy / Math.sqrt(sxx * syy);
}

export function spearman(x, y) { return pearson(ranks(x), ranks(y)); }

/** 決定的な RNG(mulberry32)。seed を固定するので再現する。 */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Spearman の両側 p を並べ替え検定で出す(|rho| 基準)。
 * ★決定的(seed 固定)。★p は (超えた数 + 1)/(反復 + 1) で 0 にならない。
 */
export function permP(x, y, iters = 20000, seed = 20260907) {
  const rx = ranks(x), ry = ranks(y);
  const r0 = Math.abs(pearson(rx, ry));
  if (!Number.isFinite(r0)) return NaN;
  const rnd = mulberry32(seed);
  const arr = ry.slice();
  let ge = 0;
  for (let it = 0; it < iters; it++) {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(rnd() * (i + 1));
      const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
    }
    if (Math.abs(pearson(rx, arr)) >= r0 - 1e-12) ge++;
  }
  return (ge + 1) / (iters + 1);
}

/** Holm-Bonferroni の調整後 p(単調化して 1 で頭打ち)。 */
export function holm(ps) {
  const m = ps.length;
  const ord = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]);
  const adj = new Array(m);
  let run = 0;
  for (let k = 0; k < m; k++) {
    const v = (m - k) * ord[k][0];
    run = Math.max(run, v);
    adj[ord[k][1]] = Math.min(1, run);
  }
  return adj;
}

/** 偏 Spearman(順位に対する 1 変量の偏相関)。 */
export function partial(rxy, rxz, ryz) {
  const d = Math.sqrt((1 - rxz * rxz) * (1 - ryz * ryz));
  if (d === 0) return NaN;
  return (rxy - rxz * ryz) / d;
}

/** 一元配置の ICC(1) と設計効果。groups: 値の配列の配列。 */
export function icc1(groups) {
  const gs = groups.filter(g => g.length > 0);
  const k = gs.length;
  const N = gs.reduce((s, g) => s + g.length, 0);
  if (k < 2 || N <= k) return { icc: NaN, deff: NaN, nEff: N, k, N };
  const grand = gs.flat().reduce((s, v) => s + v, 0) / N;
  let msb = 0, msw = 0;
  for (const g of gs) {
    const m = g.reduce((s, v) => s + v, 0) / g.length;
    msb += g.length * (m - grand) ** 2;
    for (const v of g) msw += (v - m) ** 2;
  }
  msb /= (k - 1);
  msw /= (N - k);
  const n0 = (N - gs.reduce((s, g) => s + g.length ** 2, 0) / N) / (k - 1);
  const icc = (msb - msw) / (msb + (n0 - 1) * msw);
  const mbar = N / k;
  const deff = 1 + (mbar - 1) * Math.max(0, icc);
  return { icc, deff, nEff: N / deff, k, N, n0 };
}

// ════════════════════════════════════════════════════════════════════
// 3. 読み出し
// ════════════════════════════════════════════════════════════════════

export function projectRoots(slugPrefix = 'D--Math-ABC3') {
  const base = path.join(os.homedir(), '.claude', 'projects');
  if (!fs.existsSync(base)) return [];
  return fs.readdirSync(base)
    .filter(d => d === slugPrefix || d.startsWith(slugPrefix + '-'))
    .map(d => path.join(base, d));
}

/** 通知 + 子 agent の meta を突き合わせて 1 件 = 1 agent にする。 */
export function collect(opts = {}) {
  const roots = opts.roots ?? projectRoots();
  const byToolUse = new Map();   // toolUseId -> record
  const metaByToolUse = new Map();
  const subagentFile = new Map();

  for (const root of roots) {
    if (!fs.existsSync(root)) continue;
    // (a) 子 agent の meta / 本文
    for (const sess of fs.readdirSync(root)) {
      const dir = path.join(root, sess, 'subagents');
      if (!fs.existsSync(dir)) continue;
      for (const f of fs.readdirSync(dir)) {
        if (!f.endsWith('.meta.json')) continue;
        let m; try { m = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { continue; }
        if (!m.toolUseId) continue;
        metaByToolUse.set(m.toolUseId, m);
        const jf = path.join(dir, f.replace(/\.meta\.json$/, '.jsonl'));
        if (fs.existsSync(jf)) subagentFile.set(m.toolUseId, jf);
      }
    }
    // (b) 親セッションの通知
    for (const f of fs.readdirSync(root)) {
      if (!f.endsWith('.jsonl')) continue;
      const text = fs.readFileSync(path.join(root, f), 'utf8');
      for (const line of text.split('\n')) {
        if (!line.includes('duration_ms')) continue;
        let o; try { o = JSON.parse(line); } catch { continue; }
        const content = typeof o.content === 'string' ? o.content : null;
        if (!content) continue;
        const p = parseNotification(content);
        if (!p) continue;
        const prev = byToolUse.get(p.toolUseId);
        // 同じ agent が 2 回以上通知しうる(resume)。★最後の(= 最大の)duration を採る。
        if (!prev || p.durationMs > prev.durationMs) {
          byToolUse.set(p.toolUseId, { ...p, ts: o.timestamp ?? null });
        }
      }
    }
  }

  const out = [];
  for (const [tu, rec] of byToolUse) {
    const m = metaByToolUse.get(tu) ?? {};
    out.push({
      ...rec,
      agentType: m.agentType ?? '?',
      description: m.description ?? rec.name,
      spawnDepth: m.spawnDepth ?? null,
      subagentFile: subagentFile.get(tu) ?? null,
      day: (rec.ts ?? '').slice(0, 10),
    });
  }
  out.sort((a, b) => (a.ts ?? '') < (b.ts ?? '') ? -1 : 1);
  return out;
}

/** 子 agent の本文から共変量を読む(重いので必要なときだけ)。 */
export function covariates(rec, repoRoot, opts = {}) {
  const zero = { leanChecks: 0, checkFails: 0, writes: 0, edits: 0, bashes: 0,
                 briefChars: 0, est: null, estMid: NaN, files: [], lines: NaN, cores: NaN, file: null,
                 linesWritten: NaN, coresWritten: NaN };  // ★linesWritten は「最後の Write の本文」であって完了時の行数ではない(その後の Edit を含まない)
  if (!rec.subagentFile || !fs.existsSync(rec.subagentFile)) return zero;
  const text = opts.text ?? fs.readFileSync(rec.subagentFile, 'utf8');
  return covariatesFromText(text, repoRoot, opts);
}

/**
 * ★本文(jsonl の中身)から共変量を作る純関数。★`covariates` はこれに委譲するだけ。
 * ★分けた理由(メタ第 30 回): ★**分母を突然変異させて族が動くかを試す**ため。
 *   ファイルを書き換えずに本文だけ差し替えられないと、`--denominator` の実証ができない。
 *
 * opts:
 *   fileSource : 'touched'(既定・現行) … `.lean` の path が tool_use の入力に現れた本すべて
 *                'written'             … ★その agent が Write|Edit した本だけ
 *   failRule   : 'text'(既定・現行)    … 結果の字面に /error|✗|failed/i があれば失敗
 *                'header'              … ★lean_check 自身の見出し「エラー N 件」で数える
 * ★既定は現行と 1 ビットも変わらないこと(selftest が見ている)。
 */
export function covariatesFromText(text, repoRoot, opts = {}) {
  const fileSource = opts.fileSource ?? 'touched';
  const failRule = opts.failRule ?? 'text';
  const zero = { leanChecks: 0, checkFails: 0, writes: 0, edits: 0, bashes: 0,
                 briefChars: 0, est: null, estMid: NaN, files: [], lines: NaN, cores: NaN, file: null,
                 linesWritten: NaN, coresWritten: NaN };
  const lines = text.split('\n');
  const c = { ...zero, files: [] };
  const fileHits = new Map();
  const written = new Map();
  let pendingCheck = new Set();
  for (const line of lines) {
    if (!line) continue;
    let o; try { o = JSON.parse(line); } catch { continue; }
    const msg = o.message;
    if (!msg) continue;
    const content = Array.isArray(msg.content) ? msg.content : [];
    for (const b of content) {
      if (b.type === 'tool_use') {
        const n = b.name || '';
        if (n === 'Write') c.writes++;
        else if (n === 'Edit' || n === 'MultiEdit') c.edits++;
        else if (n === 'Bash') c.bashes++;
        if (/lean_check$/.test(n)) { c.leanChecks++; pendingCheck.add(b.id); }
        const fp = b.input && (b.input.file_path || b.input.filePath);
        if (typeof fp === 'string' && /\.lean$/i.test(fp) && /Found[\\/]/.test(fp)) {
          const isW = (n === 'Write' || n === 'Edit' || n === 'MultiEdit');
          if (fileSource !== 'written' || isW) fileHits.set(fp, (fileHits.get(fp) || 0) + 1);
          // ★完了時点の行数。★ファイルは後から別の agent に書き換えられるので、
          //   「いま wc -l した行数」ではなく **その agent が書いた最後の本文** を採る。
          if (n === 'Write' && typeof b.input.content === 'string') {
            written.set(fp, b.input.content);
          }
        }
      } else if (b.type === 'tool_result') {
        if (pendingCheck.has(b.tool_use_id)) {
          pendingCheck.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          const bad = failRule === 'header' ? /エラー\s*\d+\s*件/.test(s) : /error|✗|failed/i.test(s);
          if (bad) c.checkFails++;
        }
      }
    }
    if (c.briefChars === 0 && msg.role === 'user' && typeof msg.content === 'string') {
      c.briefChars = msg.content.length;
      c.est = parseEstimate(msg.content);
      c.estMid = c.est ? (c.est.lo + c.est.hi) / 2 : NaN;
    }
  }
  c.files = [...fileHits.entries()].sort((a, b) => b[1] - a[1]).map(e => e[0]);
  if (c.files.length) {
    c.file = c.files[0];
    const w = written.get(c.file);
    if (typeof w === 'string') {
      c.linesWritten = w.split('\n').length - (w.endsWith('\n') ? 1 : 0);
      c.coresWritten = countAbstractCores(w);
    }
    const abs = path.isAbsolute(c.file) ? c.file : path.join(repoRoot, c.file);
    const cand = [abs, path.join(repoRoot, c.file.replace(/^.*?lean[\\/]/, 'lean/'))];
    for (const p2 of cand) {
      if (fs.existsSync(p2)) {
        const src = fs.readFileSync(p2, 'utf8');
        c.lines = src.split('\n').length - (src.endsWith('\n') ? 1 : 0);
        c.cores = countAbstractCores(src);
        break;
      }
    }
  }
  return c;
}

// ════════════════════════════════════════════════════════════════════
// 4. 印字
// ════════════════════════════════════════════════════════════════════
const min1 = (ms) => (ms / 60000).toFixed(1);
const pad = (s, n) => String(s).padStart(n);
const padr = (s, n) => { s = String(s); return s + ' '.repeat(Math.max(0, n - width(s))); };
function width(s) { let w = 0; for (const ch of String(s)) w += /[　-鿿＀-￯]/.test(ch) ? 2 : 1; return w; }

function describe(vals) {
  return {
    n: vals.length,
    min: Math.min(...vals), q1: quantile(vals, 0.25), med: quantile(vals, 0.5),
    q3: quantile(vals, 0.75), max: Math.max(...vals),
    mean: vals.reduce((s, v) => s + v, 0) / vals.length,
    sum: vals.reduce((s, v) => s + v, 0),
  };
}

function cmdStats(recs) {
  const d = describe(recs.map(r => r.durationMs));
  console.log(`## 1 ノードあたりの時間の分布 —— n = ${d.n}`);
  console.log('');
  console.log('  統計量        ミリ秒        分');
  const row = (k, v) => console.log(`  ${padr(k, 12)}${pad(Math.round(v), 10)}   ${pad(min1(v), 7)}`);
  row('最小', d.min); row('第1四分位', d.q1); row('中央値', d.med);
  row('第3四分位', d.q3); row('最大', d.max); row('平均', d.mean);
  console.log(`  ${padr('四分位範囲', 12)}${pad(Math.round(d.q3 - d.q1), 10)}   ${pad(min1(d.q3 - d.q1), 7)}`);
  console.log(`  ${padr('合計', 12)}${pad(Math.round(d.sum), 10)}   ${pad(min1(d.sum), 7)}`);
  console.log('');
  // 実時間の窓(重なりを含む) —— 直列/並行の差を出す
  const ts = recs.map(r => r.ts).filter(Boolean).sort();
  if (ts.length >= 2) {
    const span = new Date(ts[ts.length - 1]) - new Date(ts[0]);
    console.log(`  最初の完了 ${ts[0]} / 最後の完了 ${ts[ts.length - 1]}`);
    console.log(`  ★完了時刻の幅 ${min1(span)} 分 に対し agent 時間の合計は ${min1(d.sum)} 分`);
    console.log(`    ⇒ 重なり率 = 合計 / 幅 = ${(d.sum / span).toFixed(2)} 倍(1.0 なら直列)`);
    console.log(`    ⇒ ★「幅 / 件数」= ${min1(span / d.n)} 分/ノード(★これが 37.2 分と同じ定義の量)`);
  }
  return d;
}

function fmtP(p) { return Number.isFinite(p) ? p.toFixed(4) : '  ——  '; }

function cmdExplain(rows) {
  console.log('## ★何が時間を説明するか(★向きは断定しない)');
  console.log('');
  console.log('  ★比較の族は**データを見る前に**コードへ固定してある(FAMILY)。');
  const usable = [], dropped = [];
  for (const f of FAMILY) {
    const sub = usableRows(rows, f.key);
    if (sub.length >= MIN_N && new Set(sub.map(r => r[f.key])).size > 1) usable.push({ f, sub });
    else dropped.push({ f, sub });
  }
  console.log(`  検定できるのは ${usable.length} 本 / 族 ${FAMILY.length} 本 ⇒ Holm で補正する(m = ${usable.length})`);
  for (const { f, sub } of dropped) {
    const miss = rows.length - sub.length;
    const why = sub.length < MIN_N ? `件数不足 n=${sub.length} < ${MIN_N}` : '全件同値';
    console.log(`  ★測れなかった: ${f.label}(欠測 ${miss} / ${rows.length}、${why})`);
  }
  console.log('');
  const raw = usable.map(({ f, sub }) => {
    const yy = sub.map(r => r.durationMs);
    const xx = sub.map(r => r[f.key]);
    return { f, n: sub.length, rho: spearman(xx, yy), p: permP(xx, yy) };
  });
  const adj = holm(raw.map(r => r.p));
  console.log('  説明変数                                        n   Spearman ρ    並べ替え p   Holm 後   判定');
  raw.forEach((r, i) => {
    const verdict = adj[i] < 0.05 ? '★言える(補正後)' : '言えない';
    console.log(`  ${padr(r.f.label, 42)} ${pad(r.n, 4)} ${pad(r.rho.toFixed(3), 8)}   ${pad(fmtP(r.p), 10)}  ${pad(fmtP(adj[i]), 8)}   ${verdict}`);
  });
  console.log('');
  console.log('  ★n が欄ごとに違うのは **欄ごとの complete-case** だから(欠測の理由が欄で違う)。');
  console.log('');
  console.log('  ★「言えない」は「差が無い」ではない。★ρ の符号を根拠に向きを主張しないこと。');
  console.log('');
  // 交絡
  console.log('## ★交絡(説明変数どうしの Spearman ρ)');
  console.log('');
  const keys = usable.map(u => u.f.key);
  const pairRho = (a, b) => {
    const s2 = rows.filter(r => Number.isFinite(r[a]) && Number.isFinite(r[b]));
    return spearman(s2.map(r => r[a]), s2.map(r => r[b]));
  };
  console.log('  ' + padr('', 12) + keys.map(k => pad(k.slice(0, 10), 11)).join(''));
  for (const a of keys) {
    let line = '  ' + padr(a.slice(0, 12), 12);
    for (const b of keys) {
      const v = a === b ? 1 : pairRho(a, b);
      line += pad(Number.isFinite(v) ? v.toFixed(3) : '——', 11);
    }
    console.log(line);
  }
  console.log('');
  // 偏相関(族の外。断定しない)
  if (keys.includes('toolUses') && keys.includes('tokens')) {
    const s3 = rows.filter(r => Number.isFinite(r.toolUses) && Number.isFinite(r.tokens));
    const y = s3.map(r => r.durationMs);
    const t = s3.map(r => r.toolUses), k = s3.map(r => r.tokens);
    const rtk = spearman(t, k), ryt = spearman(y, t), ryk = spearman(y, k);
    console.log(`  ★tool_uses と tokens の ρ = ${rtk.toFixed(3)}`);
    console.log(`    duration~tool_uses を tokens で偏らせる  ρ = ${partial(ryt, ryk, rtk).toFixed(3)}`);
    console.log(`    duration~tokens を tool_uses で偏らせる  ρ = ${partial(ryk, ryt, rtk).toFixed(3)}`);
    console.log('    ★これは族の外(検定していない)。★交絡の大きさを見るためだけに出す。');
  }
  console.log('');
  return { raw, adj };
}

function bundleKey(desc) {
  // 「B2+B3 …」「段1c …」「Λ9 …」等、先頭の記号列を束の鍵にする(枝番は落とす)
  const head = String(desc).trim().split(/[\s—]/)[0].split('+')[0];
  return head.replace(/[0-9]+[a-z]*$/, '') || head;
}

function cmdIcc(rows) {
  const g = new Map();
  for (const r of rows) {
    const k = bundleKey(r.description);
    if (!g.has(k)) g.set(k, []);
    g.get(k).push(r.durationMs);
  }
  const groups = [...g.values()];
  const s = icc1(groups);
  console.log('## ★標本は独立か —— 束(持ち場の頭)ごとのばらつき');
  console.log('');
  console.log('  ' + [...g.entries()].sort((a, b) => b[1].length - a[1].length)
    .map(([k, v]) => `${k}:${v.length}`).join(' '));
  console.log(`  束 ${s.k} / 件数 ${s.N} / ICC = ${Number.isFinite(s.icc) ? s.icc.toFixed(3) : '——'}` +
              ` / DEFF = ${Number.isFinite(s.deff) ? s.deff.toFixed(2) : '——'}` +
              ` ⇒ 有効件数 ≒ ${Number.isFinite(s.nEff) ? s.nEff.toFixed(0) : '——'}`);
  console.log('  ★束が 1 件ずつなら ICC は意味を持たない(そのときは「粗い」と読むこと)。');
  console.log('');
  return s;
}

function cmdEstimate(rows) {
  const withEst = rows.filter(r => r.est && Number.isFinite(r.lines));
  console.log('## ★見積(brief の「見積 N–M 行」)と実測(wc -l)');
  console.log('');
  if (withEst.length === 0) { console.log('  ★見積を書いた brief が 0 件。測れない。'); return null; }
  console.log('  持ち場                                     見積      実測   実測/見積中点   範囲内');
  let inRange = 0;
  for (const r of withEst) {
    const mid = (r.est.lo + r.est.hi) / 2;
    const ok = r.lines >= r.est.lo && r.lines <= r.est.hi;
    if (ok) inRange++;
    console.log(`  ${padr(String(r.description).slice(0, 34), 36)} ${pad(r.est.lo + '-' + r.est.hi, 9)} ${pad(r.lines, 7)}   ${pad((r.lines / mid).toFixed(2), 10)}    ${ok ? '○' : '×'}`);
  }
  const ratio = withEst.map(r => r.lines / ((r.est.lo + r.est.hi) / 2));
  console.log('');
  console.log(`  n = ${withEst.length} / 範囲内 ${inRange} 件(${(100 * inRange / withEst.length).toFixed(0)}%)`);
  console.log(`  実測/見積中点: 中央値 ${quantile(ratio, 0.5).toFixed(2)} / 四分位 ${quantile(ratio, 0.25).toFixed(2)}–${quantile(ratio, 0.75).toFixed(2)} / 範囲 ${Math.min(...ratio).toFixed(2)}–${Math.max(...ratio).toFixed(2)}`);
  console.log('  ★これは「行数の見積」であって「かかった時間の見積」ではない。');
  return { n: withEst.length, inRange, ratio };
}

function cmdCost(rows, paths) {
  const { costs, missing } = readCosts(paths);
  console.log('## ★自己申告(COST[<持ち場>]: 安|並|高)と実測');
  console.log('');
  for (const m of missing) console.log(`  (記録が無い: ${m})`);
  console.log(`  読んだ記録: ${paths.filter(p => !missing.includes(p)).join(' / ') || '(無し)'}`);
  console.log(`  申告 ${costs.size} 件 / 実測できる持ち場 ${rows.length} 件`);
  if (costs.size === 0) {
    console.log('');
    console.log('  ★申告が 0 件。**分母が立たないので何も言わない。**');
    console.log('  ★書式:  COST[<持ち場>]: <安|並|高>  — <一言>   (ResearchPaper/decisions-pending.md)');
    console.log('  ★「安い/高い」の語を本文から拾う数え方は M96 で**測れない**と分かっている');
    console.log('    (695 節中 21 節が語を含むのに、持ち場に紐づいたのは 1 件だけだった)。');
    return { n: 0, joined: 0, costs };
  }
  const joined = [];
  const unjoined = [];
  for (const [k, c] of costs) {
    const hits = rows.filter(r => matchesKey(r.description, c.slot ?? k));
    if (hits.length === 0) { unjoined.push([k, c]); continue; }
    for (const r of hits) joined.push({ k, c, r });
  }
  console.log('');
  console.log('  持ち場                              申告   所要(分)   実測行   見積中点  実測/見積');
  for (const j of joined) {
    const mid = j.r.est ? (j.r.est.lo + j.r.est.hi) / 2 : NaN;
    console.log(`  ${padr(String(j.r.description).slice(0, 30), 32)} ${padr(j.c.level, 5)} ${pad(min1(j.r.durationMs), 8)} ${pad(Number.isFinite(j.r.lines) ? j.r.lines : '——', 8)} ${pad(Number.isFinite(mid) ? mid : '——', 10)} ${pad(Number.isFinite(mid) && Number.isFinite(j.r.lines) ? (j.r.lines / mid).toFixed(2) : '——', 9)}`);
  }
  console.log('');
  const tally = COST_LEVELS.map(l => `${l}:${joined.filter(j => j.c.level === l).length}`).join(' / ');
  console.log(`  突き合わせ ${joined.length} 件(${tally})/ 実測に当たらなかった申告 ${unjoined.length} 件`);
  for (const [k, c] of unjoined) console.log(`    ・${k}(${c.level}) —— ${c.file}:${c.line}`);
  if (unjoined.length) {
    console.log('    ⇒ ★鍵が agent の呼び名と違うだけなら `| 持ち場=<呼び名>` を 1 つ足せば繋がる。');
  }
  const usable = joined.filter(j => j.r.est && Number.isFinite(j.r.lines));
  console.log('');
  console.log(`  ★見積が読めて実測もある組 ${usable.length} 件。`);
  const levels = usable.map(j => j.c.level);
  const tally2 = COST_LEVELS.map(l => `${l}:${levels.filter(x => x === l).length}`).join(' / ');
  if (!costFamilyReady(levels)) {
    console.log(`  ★★族に入れる条件を満たさない(内訳 ${tally2})。**判定を出さない。★向きも書かない。**`);
    console.log(`    ⇒ 事前登録した条件は 2 つ: (a) ${MIN_N} 件以上 (b) ★各水準に 3 件以上。`);
    console.log('    ⇒ ★(b) が要るのは、件数だけ足りても「全部 安」なら分散ゼロで順位相関が定義できないため。');
    console.log('    ⇒ 両方満たしたときに **1 度だけ**「申告 × 実測/見積」を族に足すこと(m が 7 → 8 になる)。');
  } else {
    const lv = usable.map(j => COST_LEVELS.indexOf(j.c.level));
    const rt = usable.map(j => j.r.lines / ((j.r.est.lo + j.r.est.hi) / 2));
    const p = permP(lv, rt);
    console.log(`  申告(安=0/並=1/高=2)と 実測/見積 の Spearman ρ = ${spearman(lv, rt).toFixed(3)} / 並べ替え p = ${fmtP(p)}`);
    console.log('  ★これは**族の外**(この道具が初めて出す欄)。★次の波で FAMILY に入れて確かめること。');
  }
  return { n: costs.size, joined: joined.length, costs };
}

// ════════════════════════════════════════════════════════════════════
// 4.5 --denominator —— ★「数えているのは実行か、書かれた字面か」を欄ごとに確かめる
// ════════════════════════════════════════════════════════════════════
/**
 * ★M142(メタ第 29 回が残した宿題)の口。
 * ★M136 が `lake build` で見つけた取り違え —— **記録の本文に書き写した命令を実行と数える** ——
 *   が、事前登録した族 FAMILY にも入っているかを、★**実データで**確かめる。
 * ★判定を差し替えるためのものではない。**並べて出すだけ**。
 */

/** 族 7 本の Holm の階段を返す純関数(★cmdExplain と同じ規則: MIN_N 以上 かつ 全件同値でない)。 */
export function familyLadder(rows) {
  const usable = [];
  for (const f of FAMILY) {
    const sub = usableRows(rows, f.key);
    if (sub.length >= MIN_N && new Set(sub.map(r => r[f.key])).size > 1) usable.push({ f, sub });
  }
  const raw = usable.map(({ f, sub }) => ({
    key: f.key, label: f.label, n: sub.length,
    rho: spearman(sub.map(r => r[f.key]), sub.map(r => r.durationMs)),
    p: permP(sub.map(r => r[f.key]), sub.map(r => r.durationMs)),
  }));
  const adj = holm(raw.map(r => r.p));
  raw.forEach((r, i) => { r.holm = adj[i]; r.say = adj[i] < 0.05; });
  return raw;
}

/** ★族の欄ごとの「出所」。★コードを読んで書いた表ではなく、下の突然変異が裏を取る。 */
export const PROVENANCE = [
  { key: 'lines',      from: '実行(tool_use の入力 file_path)+ 円盤の中身', note: '★どの本を採るかは「path が現れた回数」なので Read も混ざる' },
  { key: 'toolUses',   from: '実行(完了通知 <tool_uses>)',                  note: '系が数えた値。字面ではない' },
  { key: 'tokens',     from: '実行(完了通知 <subagent_tokens>)',            note: '系が数えた値。字面ではない' },
  { key: 'cores',      from: '実行(同上)+ 円盤の中身',                      note: '★lines と同じ本を見る' },
  { key: 'leanChecks', from: '実行(tool_use の name)',                      note: '本文には作れない' },
  { key: 'checkFails', from: '★字面(tool_result の中身 /error|✗|failed/i)', note: '★道具自身の見出しと突き合わせる' },
  { key: 'estMid',     from: '字面(brief の本文。★意図的)',                 note: '見積は書かれた字面そのもの' },
];

/** ★突然変異: Bash の命令文を「罠語だらけ」に置換した本文を作る。★族が動けば字面を読んでいる証拠。 */
export function trapMutate(text, trap) {
  const out = [];
  let n = 0;
  for (const line of text.split('\n')) {
    if (!line || !line.includes('"Bash"')) { out.push(line); continue; }
    let o; try { o = JSON.parse(line); } catch { out.push(line); continue; }
    const cs = Array.isArray(o.message?.content) ? o.message.content : [];
    let hit = false;
    for (const b of cs) if (b.type === 'tool_use' && b.name === 'Bash' && b.input && typeof b.input.command === 'string') { b.input.command = trap; hit = true; }
    if (!hit) { out.push(line); continue; }
    n++; out.push(JSON.stringify(o));
  }
  return { text: out.join('\n'), n };
}

export const TRAP_CMD = "cat > log.md <<'EOF'\nlake build ABC3 && node tools/check.mjs && node tools/graph.mjs\nlean/ABC3/Found/ZZZPhantom.lean error ✗ failed lean_check\nEOF";

function cmdDenominator(recs, repoRoot) {
  const KEYS = ['lines', 'cores', 'leanChecks', 'checkFails', 'estMid', 'writes', 'edits', 'bashes'];
  console.log('## ★1. 族の 7 欄は「実行」を数えているか「書かれた字面」を数えているか');
  console.log('');
  console.log('  欄            出所                                            覚え書き');
  for (const p of PROVENANCE) console.log(`  ${padr(p.key, 13)} ${padr(p.from, 46)} ${p.note}`);
  console.log('');

  // -- 突然変異: Bash の命令文をすべて罠語に置き換える
  let mutated = 0, moved = 0;
  const movedRows = [];
  for (const r of recs) {
    if (!r.subagentFile || !fs.existsSync(r.subagentFile)) continue;
    const src = fs.readFileSync(r.subagentFile, 'utf8');
    const m = trapMutate(src, TRAP_CMD);
    mutated += m.n;
    if (!m.n) continue;
    const a = covariates(r, repoRoot);
    const b = covariates(r, repoRoot, { text: m.text });
    const d = KEYS.filter(k => String(a[k]) !== String(b[k]));
    if (d.length) { moved++; movedRows.push([String(r.description).slice(0, 32), d.join(',')]); }
  }
  console.log('## ★2. 突然変異 —— Bash の命令文を丸ごと「罠語」に置き換える');
  console.log('');
  console.log(`  置換した命令 ${mutated} 件 / ★共変量が動いた agent ${moved} / ${recs.length}`);
  for (const x of movedRows.slice(0, 10)) console.log(`    ${x[0]} : ${x[1]}`);
  console.log(moved === 0
    ? '  ⇒ ★**族の 7 欄は Bash の命令文を 1 文字も読んでいない**(M136 の heredoc 取り違えは族に入っていない)。'
    : '  ⇒ ★★命令文の字面が族に効いている。上の欄を疑うこと。');
  console.log('');

  // -- lines/cores の分母: 読んだだけの本か、自分が書いた本か
  let nT = 0, nW = 0, onlyRead = 0;
  for (const r of recs) {
    const a = covariates(r, repoRoot);
    const b = covariates(r, repoRoot, { fileSource: 'written' });
    if (a.file) nT++;
    if (b.file) nW++;
    if (a.file && a.file !== b.file) onlyRead++;
  }
  console.log('## ★3. lines / cores の分母 —— 「触った本」か「自分が書いた本」か');
  console.log('');
  console.log(`  対象が決まった agent: 触った本 ${nT} / 自分が書いた本 ${nW} / ★代表が入れ替わる ${onlyRead}`);
  console.log('  ★入れ替わる分は「読んだだけの本の行数」を仕事の大きさとして数えていた。');
  console.log('');

  // -- checkFails の字面判定は、道具自身の見出しと何件食い違うか
  let nChk = 0, byText = 0, byHead = 0, dis = 0;
  for (const r of recs) {
    if (!r.subagentFile || !fs.existsSync(r.subagentFile)) continue;
    const pend = new Set();
    for (const line of fs.readFileSync(r.subagentFile, 'utf8').split('\n')) {
      if (!line) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }
      for (const b of (Array.isArray(o.message?.content) ? o.message.content : [])) {
        if (b.type === 'tool_use') { if (/lean_check$/.test(b.name || '')) pend.add(b.id); }
        else if (b.type === 'tool_result' && pend.has(b.tool_use_id)) {
          pend.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          const t = /error|✗|failed/i.test(s), h = /エラー\s*\d+\s*件/.test(s);
          nChk++; if (t) byText++; if (h) byHead++; if (t !== h) dis++;
        }
      }
    }
  }
  console.log('## ★3b. checkFails の字面判定は道具自身の見出しと合っているか');
  console.log('');
  console.log(`  lean_check の結果 ${nChk} 件 / 字面 /error|✗|failed/i ${byText} / 見出し「エラー N 件」${byHead} / ★食い違い ${dis}`);
  console.log(dis * 200 < nChk
    ? '  ⇒ ★食い違いは 0.5% 未満。★checkFails の字面判定は**動かない**(直しても判定は変わらない。下の v2)。'
    : '  ⇒ ★★食い違いが大きい。checkFails の分母を疑うこと。');
  console.log('');

  // -- 階段
  const mk = (opts) => recs.map(r => ({ ...r, ...covariates(r, repoRoot, opts) })).filter(r => Number.isFinite(r.durationMs));
  const V = [
    ['v0 いまのまま(事前登録)', {}],
    ['v1 lines/cores を「自分が書いた本」に限る', { fileSource: 'written' }],
    ['v2 checkFails を lean_check 自身の見出しで数える', { failRule: 'header' }],
  ];
  const tab = V.map(([nm, o]) => [nm, familyLadder(mk(o))]);
  console.log('## ★4. 分母を直したときに Holm の階段はどう動くか');
  console.log('');
  for (const [nm, st] of tab) {
    console.log(`  [${nm}]`);
    console.log('   説明変数        n       ρ    並べ替え p   Holm 後   判定');
    for (const s of st) {
      console.log(`   ${padr(s.key, 13)}${pad(s.n, 5)}${pad(s.rho.toFixed(3), 9)}${pad(fmtP(s.p), 12)}${pad(fmtP(s.holm), 10)}   ${s.say ? '★言える(補正後)' : '言えない'}`);
    }
    console.log('');
  }
  const base = tab[0][1];
  console.log('  -- 動いた欄 --');
  let flips = 0;
  for (const [nm, st] of tab.slice(1)) {
    for (const s of st) {
      const b = base.find(x => x.key === s.key);
      if (!b) continue;
      if (b.say !== s.say) { flips++; console.log(`   ★${nm.slice(0, 2)} ${padr(s.key, 12)} 判定 ${b.say ? '言える' : '言えない'} → ${s.say ? '言える' : '言えない'}(Holm ${fmtP(b.holm)} → ${fmtP(s.holm)}、n ${b.n} → ${s.n})`); }
      else if (b.n !== s.n) console.log(`    ${nm.slice(0, 2)} ${padr(s.key, 12)} 判定は同じ(n ${b.n} → ${s.n}、Holm ${fmtP(b.holm)} → ${fmtP(s.holm)})`);
    }
  }
  if (!flips) console.log('   ★判定が反転した欄は無い。');
  console.log('');
  console.log('  ★★これは「2 度目の覗き」である。★v1/v2 は**データを見た後に**決めた分母なので、');
  console.log('    ★ここに出る p を事前登録の p と同じ資格で読んではいけない(メタ第 24 回 M102 と同じ)。');
  console.log('    ★言えるのは「事前登録の判定は分母の取り方に**耐えない**」ということだけ。');
}

// ════════════════════════════════════════════════════════════════════
// 4.55 --concurrency —— ★★★M166 の事前登録をそのまま実行する口(メタ第 33 回)
// ════════════════════════════════════════════════════════════════════
/** ★純関数。recs(+ 本文)から M166 の行を作る。text は id -> 本文の Map(試験で差し替えられる)。 */
export function concurrencyRows(recs, texts) {
  const iv = recs
    .filter((r) => typeof r.ts === 'string' && Number.isFinite(r.durationMs))
    .map((r) => ({ id: r.toolUseId, s: Date.parse(r.ts) - r.durationMs, e: Date.parse(r.ts), r }))
    .filter((x) => Number.isFinite(x.s) && Number.isFinite(x.e));
  const ovAll = overlapStats(iv);
  const ivImpl = iv.filter((x) => x.r.agentType === 'lean-prover');
  const ovImpl = overlapStats(ivImpl);

  const rows = iv.map((x) => {
    const t = texts.get(x.id);
    const o = t ? concurrencyOutcomesFromText(t) : null;
    const a = ovAll.get(x.id), m = ovImpl.get(x.id) ?? null;
    return {
      ...x.r, s: x.s, e: x.e,
      allAtStart: a.atStart, allMax: a.max, allMean: a.mean,
      implAtStart: m ? m.atStart : NaN, implMax: m ? m.max : NaN, implMean: m ? m.mean : NaN,
      mcpInfra: o ? o.mcpInfra : NaN,
      mcpV1: o ? o.mcpV1 : NaN,
      busy: o ? o.kind.busy : NaN,
      timeoutK: o ? o.kind.timeout : NaN,
      nostart: o ? o.kind.nostart : NaN,
      mcpCalls: o ? o.mcpCalls : 0,
      leanStarts: o ? o.leanStarts : NaN,
      lakeCalls: o ? o.lakeCalls : 0,
      lakeWaitMs: o ? o.lakeWaitMs : NaN,
      wrote: o ? o.wrote : new Set(),
      msPerTool: x.r.toolUses > 0 ? x.r.durationMs / x.r.toolUses : NaN,
    };
  });
  // (b) 重なる相手と同じ本を書いた数。★区間が重なることが条件。
  for (const A of rows) {
    let n = 0;
    for (const f of A.wrote) {
      for (const B of rows) {
        if (B === A) continue;
        if (B.e < A.s || B.s > A.e) continue;
        if (B.wrote.has(f)) { n++; break; }
      }
    }
    A.sharedFile = n;
  }
  return rows;
}

/**
 * ★★★事後(探索)—— 「REPL は処理中」が**起きた瞬間**に、MCP を使う agent が何本走っていたか。
 * ★これは事前登録に無い。★メタ第 33 回が (a) の検出器を直したあとで思いついた測り方である。
 *   ⇒ ★**確証ではない。**★ただし「取り合いか否か」を直接見るのはこの量だけである。
 * ★母集団を lean-prover に絞らない: ★**REPL は親セッションとも共有される**(#236)。
 * 返り値 [{ t, id, agents, mcpAgents }]  ★mcpAgents は「その瞬間 MCP を 1 度でも使う agent の本数」。
 */
export function busyEvents(rows, texts) {
  const usesMcp = new Map();
  for (const r of rows) usesMcp.set(r.toolUseId, (r.mcpCalls || 0) > 0);
  const out = [];
  for (const r of rows) {
    const t = texts.get(r.toolUseId);
    if (!t) continue;
    const pend = new Set();
    for (const ln of t.split('\n')) {
      if (!ln) continue;
      let o; try { o = JSON.parse(ln); } catch { continue; }
      const cs = Array.isArray(o.message?.content) ? o.message.content : [];
      for (const b of cs) {
        if (b.type === 'tool_use' && String(b.name || '').startsWith('mcp__abc3-lean__')) pend.add(b.id);
        else if (b.type === 'tool_result' && pend.has(b.tool_use_id)) {
          pend.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          if (!/REPL は処理中/.test(s)) continue;
          const at = o.timestamp ? Date.parse(o.timestamp) : NaN;
          if (!Number.isFinite(at)) continue;
          let agents = 0, mcpAgents = 0;
          for (const q of rows) {
            if (q.s <= at && at <= q.e) { agents++; if (usesMcp.get(q.toolUseId)) mcpAgents++; }
          }
          out.push({ t: new Date(at).toISOString(), id: r.toolUseId, desc: r.description, agents, mcpAgents });
        }
      }
    }
  }
  return out.sort((a, b) => (a.t < b.t ? -1 : 1));
}

/**
 * ★★★★★事後(探索)—— **無音の環境すり替わり**を構造で見つける。
 *
 * ★D27 の訂正(decisions-pending 第 1077)が名指しした、いちばん重い壊れ方:
 *   > `lean_start(...)` は **10.2 秒**で「成功」と返ったが、`lean_status` の imports が
 *   > **もう 1 体の agent のもの**に差し替わっていた。
 * ★★これは **error の字面をひとつも出さない**。⇒ (a) の検出器では**原理的に捕まらない**。
 *   ⇒ ★代わりに **「自分が頼んだ imports」と「道具が報告した imports」の食い違い**で見る。
 *
 * ★判定(★保守的に倒す。誤報を出さないほうを選ぶ):
 *   - その agent 自身の `lean_start` の入力 `imports` の集合を**すべて**覚える。
 *   - `lean_start` / `lean_status` の出力の `imports: A, B, C` を集合にする。
 *   - ★**自分が頼んだどの集合とも一致しない**なら食い違い 1 件。
 *   - ★まだ 1 度も `lean_start` を呼んでいない agent の `lean_status` は**数えない**
 *     (誰の環境かを問う資格が無い。★ここを数えると誤報が出る)。
 * 返り値 [{ t, id, desc, got:[], want:[[...]] }]
 */
export function envMismatchEvents(rows, texts) {
  const out = [];
  const setKey = (a) => [...new Set(a)].map((s) => s.trim()).filter(Boolean).sort().join('|');

  // ── ★1 周目: 「誰がどの imports を頼んだか」の全体表を作る(★相手を名指しするため)
  const askedBy = new Map();   // setKey -> [{ id, desc, s, e }]
  for (const r of rows) {
    const text = texts.get(r.toolUseId);
    if (!text) continue;
    for (const ln of text.split('\n')) {
      if (!ln) continue;
      let o; try { o = JSON.parse(ln); } catch { continue; }
      const cs = Array.isArray(o.message?.content) ? o.message.content : [];
      for (const b of cs) {
        if (b.type === 'tool_use' && /lean_start$/.test(String(b.name || ''))) {
          const k = setKey(Array.isArray(b.input?.imports) ? b.input.imports : []);
          if (!k) continue;
          if (!askedBy.has(k)) askedBy.set(k, []);
          askedBy.get(k).push({ id: r.toolUseId, desc: r.description, s: r.s, e: r.e });
        }
      }
    }
  }

  for (const r of rows) {
    const text = texts.get(r.toolUseId);
    if (!text) continue;
    const asked = new Set();
    const askedRaw = [];
    const pend = new Map();
    for (const ln of text.split('\n')) {
      if (!ln) continue;
      let o; try { o = JSON.parse(ln); } catch { continue; }
      const cs = Array.isArray(o.message?.content) ? o.message.content : [];
      for (const b of cs) {
        if (b.type === 'tool_use' && /lean_start$/.test(String(b.name || ''))) {
          const im = Array.isArray(b.input?.imports) ? b.input.imports : [];
          asked.add(setKey(im)); askedRaw.push(im);
          pend.set(b.id, 'start');
        } else if (b.type === 'tool_use' && /lean_status$/.test(String(b.name || ''))) {
          pend.set(b.id, 'status');
        } else if (b.type === 'tool_result' && pend.has(b.tool_use_id)) {
          const kind = pend.get(b.tool_use_id); pend.delete(b.tool_use_id);
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          const m = /imports:\s*([^\\\n"]+)/.exec(s);
          if (!m) continue;
          const got = m[1].split(',').map((x) => x.trim()).filter(Boolean);
          if (!got.length) continue;
          if (asked.size === 0) continue;            // ★自分で起動していないなら問わない
          if (asked.has(setKey(got))) continue;      // ★自分が頼んだ形と一致
          // ★★相手を名指しできるか: 同じ imports を頼んだ**別の** agent が、この瞬間に走っていたか。
          const at = o.timestamp ? Date.parse(o.timestamp) : NaN;
          const cands = (askedBy.get(setKey(got)) ?? []).filter((x) => x.id !== r.toolUseId);
          const culprit = Number.isFinite(at)
            ? (cands.find((x) => x.s <= at && at <= x.e) ?? null) : null;
          out.push({ t: o.timestamp ?? null, id: r.toolUseId, desc: r.description, kind,
                     got, want: askedRaw.map((a) => a.join(', ')),
                     culprit: culprit ? culprit.desc : null,
                     anyOwner: cands.length > 0 });
        }
      }
    }
  }
  return out.sort((a, b) => (String(a.t) < String(b.t) ? -1 : 1));
}

// ══════════════════════════════════════════════════════════════════════════
// ★★M173 —— 「MCP の名指し」(2026-09-08 の規約変更)は効いたか。★見張るだけの口
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★★族を増やさない。★分子は **第 33 回の `envMismatchEvents` をそのまま**使う。
 * ★足したのは **分母**(照合の回数)だけである。分母が無いと「減った」が言えない
 *   ——★件数だけを見ると「agent の本数が減っただけ」と区別できない。
 *
 * ★分母の定義: 「**自分で `lean_start` を呼んだ** agent が受け取った、
 *   `imports:` を含む返答」1 つを 1 回の照合と数える。
 *   ⇒ `envMismatchEvents` が食い違いを探している母集団と**同じ**である(定義を合わせてある)。
 */
export function envCheckEvents(rows, texts) {
  const out = [];
  for (const r of rows) {
    const text = texts.get(r.toolUseId);
    if (!text) continue;
    let started = false;
    const pend = new Map();
    for (const ln of text.split('\n')) {
      if (!ln) continue;
      let o; try { o = JSON.parse(ln); } catch { continue; }
      const cs = Array.isArray(o.message?.content) ? o.message.content : [];
      for (const b of cs) {
        if (b.type === 'tool_use' && /lean_start$/.test(String(b.name || ''))) { started = true; pend.set(b.id, 1); }
        else if (b.type === 'tool_use' && /lean_status$/.test(String(b.name || ''))) pend.set(b.id, 1);
        else if (b.type === 'tool_result' && pend.has(b.tool_use_id)) {
          pend.delete(b.tool_use_id);
          if (!started) continue;
          const s = typeof b.content === 'string' ? b.content : JSON.stringify(b.content ?? '');
          if (/imports:\s*[^\\\n"]+/.test(s)) out.push({ t: o.timestamp ?? null, id: r.toolUseId, desc: r.description });
        }
      }
    }
  }
  return out.sort((a, b) => (String(a.t) < String(b.t) ? -1 : 1));
}

/** ★規約が変わった時刻。★**機械で決めた**: 本体の `.claude/agents/lean-prover.md` の mtime
 *  (= 「あなたは MCP を使う側か、使わない側か」を書き込んだ瞬間)。★2026-09-08 03:16:13 JST。
 *  ★書き換えたら selftest が鳴る。 */
export const MCP_WATCH_CUTOFF = '2026-09-07T18:16:13.331Z';

/** ★0 件のまま基準率 p0 を片側 α で棄却するのに要る照合の回数。
 *  ★(1 − p0)^n ≤ α ⇔ n ≥ ln α / ln(1 − p0)。★p0 が 0 か 1 以上なら**判定不能**(null)。 */
export function binomNeedZero(p0, alpha = 0.05) {
  if (!(p0 > 0) || p0 >= 1) return null;
  return Math.ceil(Math.log(alpha) / Math.log(1 - p0));
}

/**
 * ★★事前登録した判定（★データを見る前に規則を書いた）:
 *   p0     = 変更**前**の実測率(食い違い / 照合)。
 *   need   = binomNeedZero(p0)。
 *   ・変更後に 1 件でも出たら → **`まだ起きている`**（★「減った」とは言わない）
 *   ・0 件のまま n ≥ need    → **`言える`**（率は p0 より低い。片側 α=0.05）
 *   ・それ以外               → **`まだ言えない`** + あと何件かを出す
 * ★★「減った」と言えるのは 3 番目ではなく 2 番目だけである。
 */
export function mcpWatchVerdict(c, alpha = 0.05) {
  const { beforeChecks, beforeMiss, afterChecks, afterMiss } = c;
  const p0 = beforeChecks > 0 ? beforeMiss / beforeChecks : null;
  const need = p0 === null ? null : binomNeedZero(p0, alpha);
  const lines = [];
  let level;
  if (p0 === null) { level = 'まだ言えない'; lines.push('★変更前の標本が無い。基準率が作れない。'); return { level, p0, need, short: null, lines }; }
  if (afterMiss > 0) {
    level = 'まだ起きている';
    lines.push(`★★変更の後にも **${afterMiss} 件**出ている(照合 ${afterChecks} 回)。★規約は破られている。`);
    return { level, p0, need, short: null, lines };
  }
  if (afterChecks >= need) {
    level = '言える';
    lines.push(`★★0 件のまま照合 ${afterChecks} 回(必要 ${need})。★率は ${(100 * p0).toFixed(1)}% より低い(片側 α=${alpha})。`);
    return { level, p0, need, short: 0, lines };
  }
  const short = need - afterChecks;
  level = 'まだ言えない';
  lines.push(`★照合 ${afterChecks} 回 / 必要 ${need} 回。★★**あと ${short} 回**。★いまは何も言えない。`);
  return { level, p0, need, short, lines };
}

/**
 * ★★純関数。★**「このままでは溜まらない」を口が自分で言う**ための量(メタ第 37 回・持ち場 4)。
 *   ★背景: `--mcp-watch` は第 34・35・36 回と **3 セッション続けて「あと 26 回」から動いていない**。
 *   ★理由は「いまの波が MCP を使わない側に振られている」ことで、★これは**観測すれば言える**のに
 *   ★口は「あと 26 回」としか言わないので、★**人が「もうすぐ溜まる」と誤読しうる**。
 *   ⇒ ★変更後の**実測の速さ**(件/日)から、★残りに要る日数を出し、★溜まっていないならそう言う。
 *   ★これは判定(検定)ではない。★**溜まり方の観測**である。
 * 引数はすべて ISO 文字列 / 数。返り値の `level` は 3 通り:
 *   'reached'  … もう必要件数に届いている(この口では速さを言わない)
 *   'accruing' … 変更後に照合が増えており、残り日数が出せる
 *   'stalled'  … ★変更後の照合が増えていない(速さ 0)⇒ ★**このままでは溜まらない**
 */
export const MCP_ACCRUAL_STALL_DAYS = 60;
/** ★変更後の照合がこれ未満なら **速さを出さない**(M174 が `perDay` でやったのと同じ守り)。 */
export const MCP_ACCRUAL_MIN_N = 5;
export function mcpAccrual({ cut, lastTs, lastCheckTs = null, afterChecks, short }) {
  const t0 = Date.parse(cut), t1 = Date.parse(lastTs ?? '');
  const days = Number.isFinite(t0) && Number.isFinite(t1) && t1 > t0 ? (t1 - t0) / 86400000 : 0;
  const tc = Date.parse(lastCheckTs ?? '');
  const droughtH = Number.isFinite(tc) && Number.isFinite(t1) && t1 > tc ? (t1 - tc) / 3600000 : null;
  const say = [];
  if (short === 0 || short === null) return { level: 'reached', days, perDay: null, needDays: null, droughtH, lines: [] };
  if (afterChecks === 0 || days <= 0) {
    say.push(`★★変更後の照合は **${afterChecks} 回 / ${days.toFixed(2)} 日**。★速さは 0 である。`);
    say.push(`★★**このままでは溜まらない。** ★「あと ${short} 回」は待てば来る数ではない。`);
    say.push('★この口が動くには「規約変更より後に自分で `lean_start` を呼ぶ agent」が要る。');
    return { level: 'stalled', days, perDay: 0, needDays: null, droughtH, lines: say };
  }
  if (afterChecks < MCP_ACCRUAL_MIN_N) {
    // ★★ここが要 —— 少ない標本から速さを外挿しない。★時間の事実だけを言う。
    say.push(`★変更後の照合は **${afterChecks} 回**(必要 ${MCP_ACCRUAL_MIN_N} 回未満)。★★**速さ(件/日)は出さない**。`);
    say.push(`  ★理由: 変更の直後の数件から外挿すると「あと ${(days / Math.max(1, afterChecks) * short).toFixed(1)} 日」のような`
      + '嘘の安心が出る(M171 / M174 が `perDay` で踏んだ形)。');
    if (droughtH !== null) say.push(`★時間の事実だけ言う: **最後の照合から ${droughtH.toFixed(1)} 時間**、新しい照合は 1 件も無い。`);
    say.push(`★★**このままでは溜まらない**と読むべきである。★「あと ${short} 回」は待てば来る数ではない。`);
    say.push('★この口が動くには「規約変更より後に自分で `lean_start` を呼ぶ agent」が要る。');
    return { level: 'stalled', days, perDay: null, needDays: null, droughtH, lines: say };
  }
  const perDay = afterChecks / days;
  const needDays = short / perDay;
  const stalled = needDays > MCP_ACCRUAL_STALL_DAYS;
  say.push(`${stalled ? '★★' : '★'}変更後 ${afterChecks} 回 / ${days.toFixed(2)} 日 = ${perDay.toFixed(2)} 件/日`
    + ` ⇒ 残り ${short} 回に **${needDays.toFixed(1)} 日**。`);
  if (stalled) say.push(`★★**このままでは溜まらない**(${MCP_ACCRUAL_STALL_DAYS} 日を超える)。★待つのではなく波の振り方を変えるしかない。`);
  return { level: stalled ? 'stalled' : 'accruing', days, perDay, needDays, droughtH, lines: say };
}

/**
 * ★★★純関数。区間の集合について「同時 k 本だった時間」を返す(時間重み)。
 * ★これが「上限が効いているか」を直接答える量である:
 *   ★上限 2 は、**同時 2 本だった時間**の間しか効かない。そこが短ければ 3 に上げても何も増えない。
 * 返り値 { total, byK: Map(k -> ms), span }
 *   total … 少なくとも 1 本走っていた時間の合計、span … 窓の端から端まで
 */
export function occupancy(iv) {
  const ev = [];
  for (const a of iv) { ev.push([a.s, +1]); ev.push([a.e, -1]); }
  ev.sort((x, y) => (x[0] - y[0]) || (x[1] - y[1]));
  const byK = new Map();
  let k = 0, prev = null, total = 0;
  for (const [t, d] of ev) {
    if (prev !== null && t > prev && k > 0) {
      byK.set(k, (byK.get(k) || 0) + (t - prev));
      total += t - prev;
    }
    k += d; prev = t;
  }
  const ts = iv.map((a) => a.s).concat(iv.map((a) => a.e)).sort((a, b) => a - b);
  return { total, byK, span: ts.length ? ts[ts.length - 1] - ts[0] : 0 };
}

/** ★純関数。族の階段(転帰 vs 露出)。★族は FAMILY_C の 4 本に固定。 */
export function concurrencyLadder(rows, expKey) {
  const usable = [];
  for (const f of FAMILY_C) {
    const sub = rows.filter((r) => Number.isFinite(r[f.key]) && Number.isFinite(r[expKey]));
    if (sub.length >= MIN_N
      && new Set(sub.map((r) => r[f.key])).size > 1
      && new Set(sub.map((r) => r[expKey])).size > 1) usable.push({ f, sub });
  }
  const raw = usable.map(({ f, sub }) => ({
    key: f.key, label: f.label, n: sub.length,
    rho: spearman(sub.map((r) => r[expKey]), sub.map((r) => r[f.key])),
    p: permP(sub.map((r) => r[expKey]), sub.map((r) => r[f.key])),
  }));
  const adj = holm(raw.map((r) => r.p));
  raw.forEach((r, i) => { r.holm = adj[i]; r.say = adj[i] < 0.05; });
  return { rows: raw, skipped: FAMILY_C.filter((f) => !usable.some((u) => u.f.key === f.key)) };
}

/** ★★M173 —— 「MCP の名指し」の前後を見る口。★見張るだけ。★階段は計算しない。 */
function cmdMcpWatch(recs) {
  const texts = new Map();
  for (const r of recs) {
    if (r.subagentFile && fs.existsSync(r.subagentFile)) texts.set(r.toolUseId, fs.readFileSync(r.subagentFile, 'utf8'));
  }
  const rows = concurrencyRows(recs, texts);
  const checks = envCheckEvents(rows, texts);
  const miss = envMismatchEvents(rows, texts);
  const cut = MCP_WATCH_CUTOFF;
  const c = {
    beforeChecks: checks.filter((e) => String(e.t) <= cut).length,
    beforeMiss: miss.filter((e) => String(e.t) <= cut).length,
    afterChecks: checks.filter((e) => String(e.t) > cut).length,
    afterMiss: miss.filter((e) => String(e.t) > cut).length,
  };
  const v = mcpWatchVerdict(c);
  console.log('## ★★M173 「MCP の名指し」は効いたか(★分子は第 33 回の検出器のまま。★足したのは分母だけ)\n');
  console.log(`  規約が変わった時刻 : ${cut}(.claude/agents/lean-prover.md の mtime)`);
  const ts = checks.map((e) => e.t).filter(Boolean).sort();
  console.log(`  ログの窓           : ${ts[0] ?? '—'} 〜 ${ts[ts.length - 1] ?? '—'}`);
  console.log(`  照合の総数         : ${checks.length}(自分で lean_start を呼んだ agent の imports 付きの返答)`);
  console.log('');
  console.log('           照合   食い違い     率');
  console.log(`  変更前 ${String(c.beforeChecks).padStart(7)}${String(c.beforeMiss).padStart(10)}   `
    + (c.beforeChecks ? (100 * c.beforeMiss / c.beforeChecks).toFixed(1) + '%' : '—'));
  console.log(`  変更後 ${String(c.afterChecks).padStart(7)}${String(c.afterMiss).padStart(10)}   `
    + (c.afterChecks ? (100 * c.afterMiss / c.afterChecks).toFixed(1) + '%' : '—'));
  console.log(`\n  ★判定: **${v.level}**`);
  for (const l of v.lines) console.log(`    ${l}`);
  // ★★溜まり方(メタ第 37 回・持ち場 4)——「あと N 回」を「もうすぐ溜まる」と読ませない
  // ★★分母は「照合の窓の終わり」ではなく **ログ全体の終わり**である。
  //   ★照合の窓で割ると、変更の 10 分後に 1 件あるだけで「135 件/日」という嘘の速さが出る
  //   (★M171 が `perDay` で踏んだのと同じ形。★実際に一度この形で書いて気づいた)。
  const logTs = recs.map((r) => r.ts).filter((x) => typeof x === 'string').sort();
  const afterTs = checks.map((e) => String(e.t)).filter((x) => x > cut).sort();
  const acc = mcpAccrual({
    cut, lastTs: logTs[logTs.length - 1] ?? null,
    lastCheckTs: afterTs[afterTs.length - 1] ?? null,
    afterChecks: c.afterChecks, short: v.short,
  });
  if (acc.lines.length) {
    console.log('\n  ★溜まり方(★判定ではなく観測):');
    for (const l of acc.lines) console.log(`    ${l}`);
  }
  console.log('\n  ★★この口は件数を見るだけで、**族を増やさない**(M166 の FAMILY_C は動かない)。');
  console.log('  ★「減った」と言えるのは **0 件のまま必要回数に届いたとき**だけである。');
}

function cmdConcurrency(recs, repoRoot) {
  const texts = new Map();
  for (const r of recs) {
    if (r.subagentFile && fs.existsSync(r.subagentFile)) {
      texts.set(r.toolUseId, fs.readFileSync(r.subagentFile, 'utf8'));
    }
  }
  const rows = concurrencyRows(recs, texts);
  const impl = rows.filter((r) => r.agentType === 'lean-prover');

  console.log('## ★★★M166 同時実行数の上限 2 は妥当か(★事前登録はコードに焼いてある)\n');
  console.log(`  母集団   : 全 agent ${rows.length} 件 / うち lean-prover **${impl.length} 件**`);
  const ts = rows.map((r) => r.s).sort((a, b) => a - b);
  if (ts.length) {
    console.log(`  ログの窓 : ${new Date(ts[0]).toISOString()} 〜 ${new Date(ts[ts.length - 1]).toISOString()}`
      + `(${((ts[ts.length - 1] - ts[0]) / 86400000).toFixed(1)} 日)`);
  }
  console.log(`  本文が読めた : ${texts.size} 件(読めない agent は (a)(b)(c) が欠測)\n`);

  // ── ★★★検出器の感度(★これを先に出す。鳴らない検出器の「0 件」は数字ではない)
  const fx = MCP_INFRA_FIXTURES.map((f) => ({ ...f, got: MCP_INFRA_RE.test(f.s), v1: MCP_INFRA_RE_V1.test(f.s) }));
  const bad = fx.filter((f) => f.got !== f.want);
  console.log('### ★★★(a) の検出器の感度(★実物の字面で毎回確かめる)\n');
  console.log('  期待  実測  v1(事前登録)  木での件数  字面');
  for (const f of fx) {
    console.log(`  ${f.want ? '鳴る' : '黙る'}  ${f.got ? '鳴る' : '黙る'}${f.got === f.want ? '  ' : '★NG'}`
      + `  ${f.v1 ? '鳴る' : '★黙る'}       ${String(f.n).padStart(4)}  ${f.s.slice(0, 46).replace(/\s+/g, ' ')}`);
  }
  const v1Tot = rows.reduce((s, r) => s + (Number.isFinite(r.mcpV1) ? r.mcpV1 : 0), 0);
  const v2Tot = rows.reduce((s, r) => s + (Number.isFinite(r.mcpInfra) ? r.mcpInfra : 0), 0);
  console.log(`\n  ★★事前登録した v1 が木で鳴った回数 : **${v1Tot} 件** / v2 : **${v2Tot} 件**`);
  if (v1Tot === 0 && v2Tot > 0) {
    console.log('  ⇒ ★★★**事前登録した (a) は感度ゼロだった。**「0 件」は「起きていない」ではない。');
    console.log('     ★以下の (a) は **事後に直した検出器**の数字であり、★**確証ではなく探索**である。');
  }
  if (bad.length) console.log(`  ★★★検出器が壊れている(${bad.length} 件が期待と違う)。★数字を信じないこと。`);
  console.log('');

  for (const [pop, set, keyStart, keyMax] of [
    ['★lean-prover(D27 が上限を掛けている母集団。★一次)', impl, 'implAtStart', 'implMax'],
    ['参考: 全 agent', rows, 'allAtStart', 'allMax'],
  ]) {
    console.log(`### ${pop}\n`);
    console.log('  水準は「自分が**始まった瞬間**に走っていた本数」(★自分を含む)');
    console.log('  水準  体数   (a)配管 /体   ★取合い /体   (b)同じ本  (c)lake中央s (d)1toolあたりs  tool_uses中央');
    const bins = new Map();
    for (const r of set) {
      const k = Number.isFinite(r[keyStart]) ? concBin(r[keyStart]) : null;
      if (k === null) continue;
      if (!bins.has(k)) bins.set(k, []);
      bins.get(k).push(r);
    }
    for (const k of [...bins.keys()].sort((a, b) => a - b)) {
      const g = bins.get(k);
      const sum = (key) => g.reduce((s, r) => s + (Number.isFinite(r[key]) ? r[key] : 0), 0);
      const med = (key) => {
        const v = g.map((r) => r[key]).filter(Number.isFinite);
        return v.length ? quantile(v, 0.5) : NaN;
      };
      const f1 = (x) => (Number.isFinite(x) ? x.toFixed(1) : '—');
      const rate = (key) => (sum(key) / g.length).toFixed(2);
      console.log(`  ${String(k === 4 ? '4+' : k).padStart(3)} ${String(g.length).padStart(5)}`
        + `${String(sum('mcpInfra')).padStart(9)} ${rate('mcpInfra').padStart(5)}`
        + `${String(sum('busy')).padStart(9)} ${rate('busy').padStart(5)}`
        + `${String(sum('sharedFile')).padStart(11)}`
        + `${f1(med('lakeWaitMs') / 1000).padStart(13)}${f1(med('msPerTool') / 1000).padStart(15)}`
        + `${f1(med('toolUses')).padStart(15)}`);
    }
    console.log('');
    console.log('  ★★(a)(b) は**合計件数**、(c)(d) は中央値。★件数の欄は「発火が 0 か否か」だけを読むこと。\n');
  }

  // ── 事前登録した階段(一次の露出だけ)
  console.log('### 階段(Spearman ρ / 並べ替え p / Holm、m = 4。★向きは断定しない)\n');
  const lad = concurrencyLadder(impl, 'implAtStart');
  if (!lad.rows.length) {
    console.log('  ★どの欄も MIN_N に届かない、または分散ゼロ。★検定しない。');
  } else {
    console.log('  欄                                        n     ρ        p      Holm   判定');
    for (const r of lad.rows) {
      console.log(`  ${padr(r.label, 40)}${pad(r.n, 4)} ${pad(r.rho.toFixed(3), 7)} `
        + `${pad(r.p.toFixed(4), 7)} ${pad(r.holm.toFixed(4), 7)}   ${r.say ? '★言える' : '言えない'}`);
    }
  }
  for (const s of lad.skipped) console.log(`  ${padr(s.label, 40)}  —— ★件数不足か分散ゼロ。検定しない`);
  console.log('');

  // ── 事前登録した判断の規則をそのまま当てる
  const hi = impl.filter((r) => Number.isFinite(r.implAtStart) && r.implAtStart >= 3);
  const hiFire = hi.reduce((s, r) => s + (r.mcpInfra || 0) + (r.sharedFile || 0), 0);
  console.log('### ★事前登録した判断の規則をそのまま当てる\n');
  console.log(`  水準 3 以上で走った lean-prover : **${hi.length} 件**(必要 ${MIN_N} 件)`);
  if (hi.length < MIN_N) {
    console.log(`  ⇒ ★★**言えない。** あと **${MIN_N - hi.length} 件**。`);
  } else {
    console.log(`  その区画の (a)+(b) の発火 : **${hiFire} 件**`);
    console.log(hiFire === 0
      ? '  ⇒ ★3 でも (a)(b) は 1 件も発火していない(★(c)(d) は中央値の表を見ること)。'
      : '  ⇒ ★★発火している。★**2 のままにすべき**。下に内訳を出す。');
    if (hiFire) {
      for (const r of hi) {
        if ((r.mcpInfra || 0) + (r.sharedFile || 0) === 0) continue;
        console.log(`     - ${r.ts} ${padr(String(r.description).slice(0, 34), 36)}`
          + ` 水準${r.implAtStart} (a)${r.mcpInfra} (b)${r.sharedFile}`);
      }
    }
  }
  console.log('');

  // ── ★★★上限が効いているか(★これが「2 を 3 にして得があるか」の本体)
  console.log('### ★★★上限は**効いて**いるか —— 同時 k 本だった**時間**\n');
  console.log('  ★上限 2 は「同時 2 本だった時間」の間しか効かない。★そこが短ければ 3 に上げても何も増えない。\n');
  for (const [pop, set] of [['★lean-prover(D27 の上限の対象)', impl], ['参考: 全 agent', rows]]) {
    const oc = occupancy(set.map((r) => ({ s: r.s, e: r.e })));
    const h = (ms) => (ms / 3600000).toFixed(1);
    console.log(`  ${pop} —— 窓 ${h(oc.span)} 時間 / 誰かが走っていた ${h(oc.total)} 時間`
      + `(窓の ${(100 * oc.total / (oc.span || 1)).toFixed(0)}%)`);
    console.log('    同時 k    時間h     走っていた時間に占める割合');
    const ks = [...oc.byK.keys()].sort((a, b) => a - b);
    for (const k of ks) {
      console.log(`    ${String(k).padStart(6)} ${h(oc.byK.get(k)).padStart(9)}`
        + `        ${(100 * oc.byK.get(k) / (oc.total || 1)).toFixed(1)}%`);
    }
    const atCap = ks.filter((k) => k >= 2).reduce((s, k) => s + oc.byK.get(k), 0);
    console.log(`    ★同時 2 本以上だった時間 : ${h(atCap)} 時間`
      + `(走っていた時間の ${(100 * atCap / (oc.total || 1)).toFixed(1)}% / 窓の ${(100 * atCap / (oc.span || 1)).toFixed(1)}%)`);
    console.log('');
  }

  // ── ★★★事後(探索): 取り合いが起きた瞬間の同時本数
  const ev = busyEvents(rows, texts);
  console.log('### ★★★事後(探索)—— 「REPL は処理中」が起きた**瞬間**に何本走っていたか\n');
  console.log('  ★事前登録に無い測り方(検出器を直した後に思いついた)。★確証ではない。');
  console.log(`  事象 **${ev.length} 件**。★母集団は全 agent(REPL は親セッションとも共有される。#236)\n`);
  if (ev.length) {
    const hist = new Map();
    for (const e of ev) hist.set(e.mcpAgents, (hist.get(e.mcpAgents) || 0) + 1);
    console.log('  その瞬間に走っていた「MCP を使う agent」の本数   事象の数');
    for (const k of [...hist.keys()].sort((a, b) => a - b)) {
      console.log(`  ${String(k).padStart(6)}                                       ${String(hist.get(k)).padStart(6)}`);
    }
    const alone = ev.filter((e) => e.mcpAgents <= 1).length;
    console.log(`\n  ★★MCP を使う agent が **1 本以下**のときに起きた事象 : **${alone} 件 / ${ev.length}**`);
    console.log(alone
      ? '  ⇒ ★★★**取り合いの相手は agent だけではない。**★agent を 1 本に絞っても消えない事象がある\n'
        + '     (親セッション / 隔離 worktree / 別プロジェクトの REPL が同じ実体を握りうる)。'
      : '  ⇒ ★すべて agent 同士の取り合いだった。');
    console.log('');
    console.log('  内訳(新しい順に 10 件):');
    for (const e of ev.slice(-10)) {
      console.log(`     ${e.t} 走行${String(e.agents).padStart(2)} うちMCP${String(e.mcpAgents).padStart(2)}`
        + `  ${String(e.desc).slice(0, 34)}`);
    }
  }
  console.log('');

  // ── ★★★★★事後(探索): 無音の環境すり替わり(D27 訂正 第 1077 の壊れ方)
  const mm = envMismatchEvents(rows, texts);
  console.log('### ★★★★★事後(探索)—— **無音**の環境すり替わり(字面のエラーが出ない壊れ方)\n');
  console.log('  ★D27 の訂正(第 1077)が名指しした壊れ方。★`lean_start` は「成功」と返る。');
  console.log('  ★判定: 自分が頼んだ imports と、道具が報告した imports が一致しない回。');
  console.log(`  ★★食い違い **${mm.length} 件**`
    + '(★自分で lean_start を呼んでいない agent の lean_status は数えない)\n');
  if (mm.length) {
    const byAgent = new Set(mm.map((e) => e.id));
    const named = mm.filter((e) => e.culprit);
    const owned = mm.filter((e) => e.anyOwner);
    console.log(`  ★のべ ${mm.length} 件 / ★agent ${byAgent.size} 体`);
    console.log(`  ★★このうち **${named.length} 件**は「同じ imports を頼んだ**別の agent が同時に走っていた**」`);
    console.log(`  ★  さらに **${owned.length} 件**は「その imports を頼んだ agent がログの中に居る」`
      + '(★時刻は重ならないが、環境が残っていた形)');
    const byDay = new Map();
    for (const e of mm) { const d = String(e.t).slice(0, 10); byDay.set(d, (byDay.get(d) || 0) + 1); }
    console.log('  ★日ごと: ' + [...byDay.entries()].sort().map(([d, n]) => `${d} ${n} 件`).join(' / '));
    console.log('  ★★D27(2026-09-07 ユーザー承認)より**後**にも出ているかを、この行で毎回見ること。\n');
    for (const e of mm) {
      console.log(`     ${e.t} [${e.kind}] ${String(e.desc).slice(0, 30)}`);
      console.log(`        頼んだ : ${e.want[e.want.length - 1] || '(空)'}`.slice(0, 118));
      console.log(`        ★報告 : ${e.got.join(', ')}`.slice(0, 118));
      console.log(`        ★★相手: ${e.culprit ? '**' + e.culprit + '**(同時に走っていた)'
        : (e.anyOwner ? '(同じ imports を頼んだ agent は居るが同時ではない)' : '—— 特定できず')}`);
    }
  } else {
    console.log('  ⇒ ★いまの木では 1 件も出ない。★★ただし「起きていない」ではなく');
    console.log('     **この検出器で見えない**可能性を残す(下の限界を読むこと)。');
  }
  console.log('');

  // ── 交絡(★都合よく隠さない)
  console.log('### ★交絡 —— 「難しい波ほど同時本数が多い」を数字で出す\n');
  const xs = impl.filter((r) => Number.isFinite(r.implAtStart));
  const co = (k) => {
    const sub = xs.filter((r) => Number.isFinite(r[k]));
    return sub.length >= MIN_N ? spearman(sub.map((r) => r.implAtStart), sub.map((r) => r[k])) : NaN;
  };
  for (const k of ['toolUses', 'tokens', 'durationMs']) {
    const v = co(k);
    console.log(`  ρ(水準, ${padr(k, 12)}) = ${Number.isFinite(v) ? v.toFixed(3) : '—(件数不足)'}`);
  }
  const groups = new Map();
  for (const r of xs) { const d = r.day || '?'; if (!groups.has(d)) groups.set(d, []); groups.get(d).push(r.implAtStart); }
  const ic = icc1([...groups.values()]);
  console.log(`  束(日ごと) ICC(1) = ${Number.isFinite(ic.icc) ? ic.icc.toFixed(3) : '—'}`
    + ` / DEFF = ${Number.isFinite(ic.deff) ? ic.deff.toFixed(2) : '—'}`
    + ` / 実効 n = ${Number.isFinite(ic.nEff) ? ic.nEff.toFixed(1) : '—'}(生 n = ${ic.N})`);
  console.log('  ★★露出そのものが持ち場の重さと絡む。★**因果は測れない**(事前登録の (1))。\n');

  // ── 素の在庫(★件数を隠さない)
  const tot = (key, set2) => set2.reduce((s, r) => s + (Number.isFinite(r[key]) ? r[key] : 0), 0);
  console.log('### 素の在庫(全 agent)\n');
  console.log(`  MCP 呼び出し ${tot('mcpCalls', rows)} 件 / うち配管の失敗 **${tot('mcpInfra', rows)} 件**`);
  console.log(`    内訳: ★取り合い(REPL は処理中) **${tot('busy', rows)} 件** / `
    + `時間切れ ${tot('timeoutK', rows)} 件 / 環境が無い ${tot('nostart', rows)} 件`);
  console.log(`  lean_start   ${tot('leanStarts', rows)} 件`);
  console.log(`  lake build   ${tot('lakeCalls', rows)} 件`);
  console.log(`  同じ本を重なって書いた延べ **${tot('sharedFile', rows)} 件**`);
}

// ════════════════════════════════════════════════════════════════════
// 4.6 --m149 —— ★事前登録(0b)をそのまま実行する口。★停止規則を道具が強制する
// ════════════════════════════════════════════════════════════════════
/** ★M160 —— 見張り。★`--record` を渡したときだけ履歴に書く(読むだけでは汚さない)。 */
function cmdM149Watch(recs, repoRoot, record) {
  const rows = freshRows(recs)
    .map((r) => ({ ...r, ...covariates(r, repoRoot, { fileSource: M149_DENOM }) }))
    .filter((r) => Number.isFinite(r.durationMs));
  const g = m149Gate(rows);
  const cur = m149Observation(recs, g.n);
  const p = path.join(repoRoot, M149_WATCH_REL);
  let hist = [];
  if (fs.existsSync(p)) {
    try { hist = JSON.parse(fs.readFileSync(p, 'utf8')).obs || []; } catch { hist = []; }
  }
  const prev = hist.length ? hist[hist.length - 1] : null;
  const v = m149WatchVerdict(prev, cur);

  console.log('## ★M160 M149 の標本の見張り(★保存期間が最大の危険。第 31 回の名指し)');
  console.log('');
  console.log(`  履歴          : ${M149_WATCH_REL}(観測 ${hist.length} 件)`);
  console.log(`  いま          : 通知 ${cur.nRecs} / T より後 ${cur.afterT} / ★使える ${cur.usable}`
    + ` / 必要 ${M149_NEED}(あと ${Math.max(0, M149_NEED - cur.usable)})`);
  console.log(`  ログの窓      : ${cur.first ?? '?'} 〜 ${cur.last ?? '?'}`);
  if (prev) {
    console.log(`  前回          : ${prev.at} —— 通知 ${prev.nRecs} / T より後 ${prev.afterT} / 使える ${prev.usable}`);
    console.log(`                  ログの窓 ${prev.first ?? '?'} 〜 ${prev.last ?? '?'}`);
  }
  console.log('');
  const mark = { first: '  ', ok: '  ', warn: '★ ', alarm: '★★' }[v.level] || '  ';
  for (const l of v.lines) console.log(`  ${mark}${l}`);
  if (v.perDay !== null) {
    console.log('');
    console.log(`  増える速さ    : ${v.perDay.toFixed(2)} 件/日`
      + (Number.isFinite(v.etaDays) ? `  ⇒ 届くのは あと ${v.etaDays.toFixed(1)} 日` : '  ⇒ ★届かない'));
  }
  console.log('');
  if (record) {
    hist.push(cur);
    fs.mkdirSync(path.dirname(p), { recursive: true });
    fs.writeFileSync(p, JSON.stringify({ v: 1, need: M149_NEED, cutoff: M149_CUTOFF, obs: hist }, null, 1));
    console.log(`  ★履歴に書いた(観測 ${hist.length} 件目)。`);
  } else {
    console.log('  ※ 読んだだけ。履歴に残すには `--record` を足す。');
  }
  console.log('  ★★この口は「使える件数」を見るだけで、**階段を計算しない**(停止規則を破らない)。');
  return v.level;
}

function cmdM149(recs, repoRoot) {
  console.log('## ★M149 の事前登録(0b)をそのまま実行する');
  console.log('');
  console.log(`  締切 T           : ${M149_CUTOFF}(★M145 が見た標本を全部外す)`);
  console.log(`  分母             : fileSource='${M149_DENOM}'(★(A) データに依らず決めた)`);
  console.log(`  一次の欄         : ${M149_PRIMARY}`);
  console.log(`  必要件数 M149_NEED: ${M149_NEED}(両側・検出力 0.80・α=0.05/7・ρ=0.30。★DEFF>1 なので下限)`);
  console.log('');
  const fresh = freshRows(recs);
  console.log(`  通知 全部 ${recs.length} 件 / ★T より後 ${fresh.length} 件`);
  const rows = fresh
    .map(r => ({ ...r, ...covariates(r, repoRoot, { fileSource: M149_DENOM }) }))
    .filter(r => Number.isFinite(r.durationMs));
  const g = m149Gate(rows);
  console.log(`  うち ${M149_PRIMARY} が取れる(= 自分が書いた .lean がある): ★${g.n} 件`);
  console.log('');
  if (!g.open) {
    console.log(`  ★★停止規則により **階段を計算しない**。あと ★${g.short} 件 必要。`);
    console.log('');
    console.log('  ★理由: 標本は時間とともに増える。届く前に覗いてから検定すると逐次の覗きになり、');
    console.log('    第一種の過誤が膨らむ。⇒ 届くまでは**件数だけ**を見る(これが正しい結果である)。');
    console.log('  ★次の人へ: この口をもう一度叩くだけでよい。締切も分母も必要件数もコードに固定してある。');
    console.log(`    node tools/agent-timing.mjs --m149`);
    console.log('  ★★`M149_CUTOFF` / `M149_NEED` を書き換えないこと(selftest が鳴る)。');
    return;
  }
  console.log('  ★★件数が届いた。事前登録どおり族 7 本の階段を 1 度だけ計算する。');
  console.log('');
  console.log('   説明変数        n       ρ    並べ替え p   Holm 後   判定');
  for (const s of g.ladder) {
    console.log(`   ${padr(s.key, 13)}${pad(s.n, 5)}${pad(s.rho.toFixed(3), 9)}${pad(fmtP(s.p), 12)}${pad(fmtP(s.holm), 10)}   ${s.say ? '★言える(補正後)' : '言えない'}`);
  }
  const pr = g.ladder.find(s => s.key === M149_PRIMARY);
  console.log('');
  console.log(pr
    ? `  ★★一次の答え: 「抽象核の本数は所要時間を説明する」は ${pr.say ? '★立つ' : '★立たない'}(Holm ${fmtP(pr.holm)}、n ${pr.n})。`
    : `  ★一次の欄 ${M149_PRIMARY} は族の階段に載らなかった(全件同値)。`);
  console.log('  ★これは**事前登録した 1 回の検定**である。★2 度目の覗きではない。');
}

// ════════════════════════════════════════════════════════════════════
// 4b. ★★M179 —— frontier の「供給」を測る口(メタ第 35 回)
// ════════════════════════════════════════════════════════════════════
//
// ★なぜ要るか(M173 = メタ第 34 回の答え)
//   `frontier.mjs` の「着手可能 N 件」は **ファイル**を数えている。実測では
//   2026-09-07 に配った lean-prover 67 本のうち **65 本(97.0%)がその一覧に載っていない**本を
//   書いていた。★本体はそれを供給量として読んでいた。
// ★★ところが M173 の測定は**使い捨て 10 本**で行われ、M1 に従って全部消された。
//   ⇒ ★**再実行できない**。★第 34 回自身がそう書いた(M179)。ここがその口である。
//
// ★★重い。既定では**歴史の木を作らない**:
//   `--supply 2026-09-07`            … いまの木の frontier と突き合わせる(約 3 秒)
//   `--supply 2026-09-07 --history`  … ★その日の朝の commit を `git worktree add --detach` して
//                                       そこで frontier / graph を立てる(★約 40 秒。既定では走らない)
// ★★`--history` を付けないと「その日の朝の一覧」ではない。★出力にそう書く。

/** ★`.lean` のパスを graph.mjs と同じ `rel`(例 `Found/PGC/X.lean`)に直す。純関数。 */
export function toRel(p) {
  if (typeof p !== 'string') return null;
  const s = p.replace(/\\/g, '/');
  const i = s.toLowerCase().lastIndexOf('lean/abc3/');
  if (i < 0) return null;
  const r = s.slice(i + 'lean/abc3/'.length);
  return /\.lean$/i.test(r) ? r : null;
}

/** ★その agent が **Write / Edit した** `.lean` の rel。純関数。
 *  ★`covariatesFromText` の `files` を使わない理由: あちらは `Found/` に絞ってあり(M149 の族の分母)、
 *    ★`Skeleton/` を書いた agent が「木に無い」に化ける。★族には触らずに別の口を作る。 */
export function writtenLeanRels(text) {
  const out = new Set();
  for (const line of String(text).split('\n')) {
    if (!line) continue;
    let o; try { o = JSON.parse(line); } catch { continue; }
    const content = Array.isArray(o.message?.content) ? o.message.content : [];
    for (const b of content) {
      if (b.type !== 'tool_use') continue;
      if (!/^(Write|Edit|MultiEdit)$/.test(b.name || '')) continue;
      const rel = toRel(b.input?.file_path ?? b.input?.filePath);
      if (rel) out.add(rel);
    }
  }
  return [...out].sort();
}

export const SUPPLY_CLASSES = ['着手可能', '着手不可', 'sorry なし', '木に無い', '書いていない'];

/** ★1 本を frontier の一覧に照らして分類する。純関数。
 *  @param {string|null} rel   その agent の主ファイル
 *  @param {Map<string,{startable:boolean}>} frontierByRel
 *  @param {Set<string>} treeRels  その時点の木にある rel */
export function supplyClass(rel, frontierByRel, treeRels) {
  if (!rel) return '書いていない';
  const f = frontierByRel.get(rel);
  if (f) return f.startable ? '着手可能' : '着手不可';
  if (treeRels.has(rel)) return 'sorry なし';
  return '木に無い';
}

/** ★`rel` が推移的に import している startable ノード。純関数。
 *  ★M173 の「52 本が `Skeleton/PGC/Section1` ただ 1 つの下流だった」はこれで再現する。
 *  ★木に無い rel は空を返す(★「辿れない」と「上流に startable が無い」を区別しない —— 危険側)。 */
export function startableUpstream(rel, nodes, startableRels) {
  const byMod = new Map(); const byRel = new Map();
  for (const n of nodes) { byMod.set(n.mod, n); byRel.set(n.rel, n); }
  const start = byRel.get(rel);
  if (!start) return [];
  const seen = new Set([start.mod]);
  const stack = [start];
  const out = new Set();
  while (stack.length) {
    const n = stack.pop();
    for (const m of (n.imports ?? [])) {
      if (seen.has(m)) continue;
      seen.add(m);
      const nx = byMod.get(m);
      if (!nx) continue;
      if (startableRels.has(nx.rel)) out.add(nx.rel);
      stack.push(nx);
    }
  }
  return [...out].sort();
}

/** ★Lean の `import ABC3.…` を読む。純関数。★行コメントの中は取らない。 */
export function parseImports(src) {
  const out = [];
  for (const line of String(src).split('\n')) {
    const m = /^\s*import\s+([A-Za-z0-9_.]+)\s*$/.exec(line);
    if (m) out.push(m[1]);
  }
  return out;
}

/** ★★その時点の木に**無かった**本を、いまの木の import で辿れるようにする合成節点。
 *  ★なぜ要るか: 2026-09-07 に配った 71 本のうち **53 本はその朝の木に存在しない**。
  *  ★合成しないと `startableUpstream` が空を返し、★M173 の「52 本が Section1 の下流」が再現しない。
 *  ★★危険側: 借りているのは**いまの** import であって、当時のものではない。★出力にそう書く。 */
export function synthNode(rel, src) {
  return { mod: `★合成 ${rel}`, rel, imports: parseImports(src) };
}

/** ★その agent に渡した brief(最初の user メッセージ)。純関数。 */
export function briefTextOf(text) {
  for (const line of String(text).split('\n')) {
    if (!line) continue;
    let o; try { o = JSON.parse(line); } catch { continue; }
    const m = o.message;
    if (m && m.role === 'user' && typeof m.content === 'string') return m.content;
  }
  return '';
}

/** ★字面に現れる `Skeleton/…/….lean` の集合。純関数。
 *  ★★これは**辺ではない**。「その持ち場がどの節点の話か」を人が書いた字面から拾うだけである。
 *  ★import の辺で辿れない理由(実測 2026-09-07): 配った 72 本の主ファイルは `Found/PGC/*` で、
 *    ★**着手可能だった `Skeleton/PGC/Section1` はそれらを import していない**(配線がまだ無い)。 */
export function mentionedSkeletons(text) {
  const out = new Set();
  for (const m of String(text).matchAll(/Skeleton[\\/][A-Za-z0-9]+[\\/][A-Za-z0-9]+\.lean/g)) {
    out.add(m[0].replace(/\\/g, '/'));
  }
  return [...out].sort();
}

/** ★分類の集計。純関数(印字と分けてある —— M91 の教訓)。 */
export function supplyTally(rows) {
  const byClass = new Map(SUPPLY_CLASSES.map((c) => [c, 0]));
  const byAnc = new Map();
  for (const r of rows) {
    byClass.set(r.klass, (byClass.get(r.klass) ?? 0) + 1);
    if (r.klass === '着手可能' || r.klass === '着手不可') continue;
    for (const a of r.upstream ?? []) byAnc.set(a, (byAnc.get(a) ?? 0) + 1);
  }
  const n = rows.length;
  const listed = (byClass.get('着手可能') ?? 0);
  const tally = (key) => {
    const m = new Map();
    for (const r of rows) for (const x of r[key] ?? []) m.set(x, (m.get(x) ?? 0) + 1);
    return [...m.entries()].sort((a, b) => b[1] - a[1]);
  };
  return {
    n, listed,
    listedPct: n ? (100 * listed) / n : 0,
    byClass,
    byAncestor: [...byAnc.entries()].sort((a, b) => b[1] - a[1]),
    byBrief: tally('briefNodes'),
    byBody: tally('bodyNodes'),
    briefHit: rows.filter((r) => (r.briefNodes ?? []).length).length,
    bodyHit: rows.filter((r) => (r.bodyNodes ?? []).length).length,
  };
}

/** ★その木の依存グラフ(`graph.mjs --json`)。 */
function graphAt(root) {
  const G = JSON.parse(execFileSync('node', [path.join(root, 'tools', 'graph.mjs'), '--json'],
    { cwd: root, encoding: 'utf8', maxBuffer: 1 << 28 }));
  return G.nodes ?? [];
}

/** ★★歴史の木を 1 本立てて frontier / graph を読む(重い)。★立てた木は必ず外す。 */
function withTreeAt(repoRoot, day, useHistory, fn) {
  if (!useHistory) return fn(repoRoot, null);
  const commit = execFileSync(
    'git', ['rev-list', '-1', `--before=${day}T00:00:00Z`, 'master'],
    { cwd: repoRoot, encoding: 'utf8' }).trim();
  if (!commit) throw new Error(`${day} より前の commit が見つからない`);
  const tmp = path.join(os.tmpdir(), `abc3-supply-${day}-${process.pid}`);
  execFileSync('git', ['worktree', 'add', '--detach', tmp, commit], { cwd: repoRoot, stdio: 'ignore' });
  try { return fn(tmp, commit); }
  finally { try { execFileSync('git', ['worktree', 'remove', '--force', tmp], { cwd: repoRoot, stdio: 'ignore' }); } catch { /* 残っても数字は出ている */ } }
}

/** ★その木の frontier(--json)。★古い commit では旗が無いことがあるので順に落とす。 */
function frontierAt(root) {
  const tries = [
    ['--json', '--all', '--limit', '0', '--no-marks'],
    ['--json', '--all', '--limit', '0'],
    ['--json', '--all'],
    ['--json'],
  ];
  let last = null;
  for (const a of tries) {
    try {
      const s2 = execFileSync('node', [path.join(root, 'tools', 'frontier.mjs'), ...a],
        { cwd: root, encoding: 'utf8', maxBuffer: 1 << 28 });
      const j = JSON.parse(s2);
      return { rows: j.frontier ?? j, flags: a };
    } catch (e) { last = e; }
  }
  throw new Error(`frontier が立たない: ${last?.message ?? '?'}`);
}

function cmdSupply(recs, repoRoot, day, useHistory) {
  const rows0 = recs.filter((r) => r.day === day && r.agentType === DIAG_TYPE);
  console.log(`## ★★M179 frontier の「供給」(${day}、${DIAG_TYPE})\n`);
  if (!rows0.length) {
    console.log(`  その日の ${DIAG_TYPE} は 0 本。★--day を付けていると二重に絞られる。`);
    return;
  }
  const t0 = Date.now();
  const out = withTreeAt(repoRoot, day, useHistory, (root, commit) => {
    const F = frontierAt(root);
    const nodes = graphAt(root);
    // ★★当時の木に**無かった**本は、当時のグラフでは 1 歩も辿れない(実測 54 本中 52 本が届かなかった)。
    //   ⇒ ★**いまのグラフ**を借りて辿る。★借りているのは辺だけで、`startable` は当時のものを使う。
    const nodesNow = root === repoRoot ? nodes : graphAt(repoRoot);
    const nowHasRel = new Set(nodesNow.map((n) => n.rel));
    const treeRels = new Set(nodes.map((n) => n.rel));
    const frontierByRel = new Map(F.rows.map((r) => [r.rel, r]));
    const startableRels = new Set(F.rows.filter((r) => r.startable).map((r) => r.rel));
    const rows = rows0.map((r) => {
      const text = r.subagentFile && fs.existsSync(r.subagentFile) ? fs.readFileSync(r.subagentFile, 'utf8') : '';
      const rels = writtenLeanRels(text);
      const rel = rels[0] ?? null;
      const klass = supplyClass(rel, frontierByRel, treeRels);
      let use = nodes, synth = false;
      if (rel && klass === '木に無い') {
        synth = true;
        use = nodesNow;
        if (!nowHasRel.has(rel)) {
          // ★いまの木にも無い(消えた / 名前が変わった)。★せめてファイルが残っていれば合成する。
          const cur = path.join(repoRoot, 'lean', 'ABC3', rel);
          if (fs.existsSync(cur)) use = [synthNode(rel, fs.readFileSync(cur, 'utf8')), ...nodesNow];
        }
      }
      const upstream = rel ? startableUpstream(rel, use, startableRels) : [];
      const briefNodes = mentionedSkeletons(briefTextOf(text));
      const bodyNodes = mentionedSkeletons(text);
      return { rel, rels, klass, upstream, synth, briefNodes, bodyNodes };
    });
    return { rows, commit, nodes: nodes.length, frontier: F.rows.length,
             startable: startableRels.size, flags: F.flags };
  });
  const T = supplyTally(out.rows);
  console.log(`  木        : ${useHistory ? `★${day} 直前の commit ${String(out.commit).slice(0, 8)}` : '★いまの作業木(その日の朝ではない)'}`);
  console.log(`  frontier  : ${out.frontier} 節点(うち着手可能 ${out.startable})/ グラフ ${out.nodes} 節点`);
  console.log(`  旗        : ${out.flags.join(' ')}`);
  console.log(`  配った    : ${T.n} 本`);
  console.log('');
  console.log('   分類            本数     割合');
  for (const c of SUPPLY_CLASSES) {
    const v = T.byClass.get(c) ?? 0;
    console.log(`   ${padr(c, 14)}${pad(v, 5)}${pad(((100 * v) / T.n).toFixed(1) + '%', 9)}`);
  }
  console.log(`\n  ★★「着手可能」に載っていたのは ${T.listed} 本(${T.listedPct.toFixed(1)}%)。`);
  if (T.byAncestor.length) {
    console.log('\n   一覧の外の本を辿ると、上流の着手可能ノードは:');
    for (const [rel, n] of T.byAncestor.slice(0, 8)) console.log(`   ${pad(n, 5)} 本 ← ${rel}`);
    console.log('   ★★1 つの着手可能ノードから何本の持ち場が切り出せているか —— これが**供給量**である。');
    const bor = out.rows.filter((r) => r.synth && r.upstream.length).length;
    if (bor) console.log(`   ★★うち ${bor} 本は**当時の木に無かった**ので、★いまのグラフの辺を借りて辿った(借り物)。`);
  }
  console.log('\n   ── ★★帰属の水路は 3 本ある。★どれも同じ数にならない(それが答えである)。');
  const chan = [
    ['① import の辺(強い)', T.byAncestor, out.rows.filter((r) => r.upstream.length).length],
    ['② brief の字面(中)', T.byBrief, T.briefHit],
    ['③ 本文の字面(弱い)', T.byBody, T.bodyHit],
  ];
  for (const [name, list, hit] of chan) {
    const top = list[0];
    console.log(`   ${padr(name, 22)} 届いた ${pad(hit, 3)} / ${T.n} 本` +
      (top ? `   最多 ${top[1]} 本 ← ${top[0]}` : '   ——'));
  }
  console.log('   ★★①が届かないのは配線がまだ無いからで、agent が遊んでいたからではない。');
  console.log(`\n  (${((Date.now() - t0) / 1000).toFixed(1)} 秒。${useHistory ? '★歴史の木を立てて外した' : '★--history を付けると その日の朝の木で測り直す'})`);
  console.log('  ★★測れないこと: agent が**主ファイル以外**も書いた場合、ここは最初の 1 本しか見ない。');
  console.log('  ★★`書いていない` は「.lean を 1 本も Write/Edit しなかった」であって、失敗とは限らない。');
}
// ════════════════════════════════════════════════════════════════════
// 5. selftest
// ════════════════════════════════════════════════════════════════════
function selftest() {
  let ok = 0, ng = 0;
  const t = (name, cond) => { if (cond) { ok++; } else { ng++; console.log('  NG ' + name); } };
  const close = (a, b, e = 1e-9) => Math.abs(a - b) <= e;

  // --- parseNotification
  const N = '<task-id>abc</task-id>\n<tool-use-id>toolu_1</tool-use-id>\n<status>completed</status>\n' +
            '<summary>Agent "B1 テスト" finished</summary>\n<result>x</result>\n' +
            '<usage><subagent_tokens>123</subagent_tokens><tool_uses>7</tool_uses><duration_ms>4500</duration_ms></usage>';
  const p = parseNotification(N);
  t('parse: duration', p && p.durationMs === 4500);
  t('parse: tool_uses', p && p.toolUses === 7);
  t('parse: tokens', p && p.tokens === 123);
  t('parse: name から Agent "…" finished を剥がす', p && p.name === 'B1 テスト');
  t('parse: toolUseId', p && p.toolUseId === 'toolu_1');
  t('parse: usage が無ければ null', parseNotification('<summary>x</summary>') === null);
  t('parse: 文字列でなければ null', parseNotification(null) === null);

  // --- parseEstimate
  t('est: 見積 **400–700 行**', JSON.stringify(parseEstimate('見積 **400–700 行**。')) === '{"lo":400,"hi":700}');
  t('est: 見積 120-250 行', JSON.stringify(parseEstimate('見積 120-250 行（Y4）')) === '{"lo":120,"hi":250}');
  t('est: 見積 約 300 行', JSON.stringify(parseEstimate('見積 約 300 行')) === '{"lo":300,"hi":300}');
  t('est: 無ければ null', parseEstimate('行数は 400 行') === null);

  // --- countAbstractCores
  const src = ['/-! ## 0. 抽象核——一般の Galois 拡大', '/-- **抽象核 1**: … -/', '/-! ## 1. 具体層', '/-! ## 2. ★抽象核 B —— x'].join('\n');
  t('cores: 節見出しだけ数える(2)', countAbstractCores(src) === 2);
  t('cores: docstring は数えない', countAbstractCores('/-- 抽象核 1 -/') === 0);

  // --- quantile(numpy type 7 と一致)
  const q = [1, 2, 3, 4];
  t('quantile: median 2.5', close(quantile(q, 0.5), 2.5));
  t('quantile: q1 1.75', close(quantile(q, 0.25), 1.75));
  t('quantile: q3 3.25', close(quantile(q, 0.75), 3.25));
  t('quantile: 単一要素', quantile([5], 0.5) === 5);
  t('quantile: 空は NaN', Number.isNaN(quantile([], 0.5)));

  // --- ranks
  t('ranks: 同順位は平均', JSON.stringify(ranks([10, 20, 20, 30])) === '[1,2.5,2.5,4]');
  t('ranks: 逆順', JSON.stringify(ranks([3, 2, 1])) === '[3,2,1]');

  // --- pearson / spearman(既知値)
  t('pearson: 完全一致 1', close(pearson([1, 2, 3], [1, 2, 3]), 1));
  t('pearson: 完全逆 -1', close(pearson([1, 2, 3], [3, 2, 1]), -1));
  // scipy.stats.pearsonr([1,2,3,4,5],[2,1,4,3,5]) = 0.8 (scipy 1.17.1 で実測)
  t('pearson: scipy 既知値 0.8', close(pearson([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]), 0.8, 1e-12));
  // scipy.stats.spearmanr([1,2,3,4,5],[5,6,7,8,7]) = 0.8207826816681233 (scipy 1.17.1 で実測)
  t('spearman: scipy 既知値 0.8207826816681233',
    close(spearman([1, 2, 3, 4, 5], [5, 6, 7, 8, 7]), 0.8207826816681233, 1e-12));
  // scipy.stats.spearmanr([1,2,3,4,5,6],[6,5,4,3,2,1]) = -1
  t('spearman: scipy 既知値 -1', close(spearman([1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1]), -1, 1e-12));
  t('spearman: 定数は NaN', Number.isNaN(spearman([1, 1, 1], [1, 2, 3])));

  // --- mulberry32 の決定性
  const r1 = mulberry32(42), r2 = mulberry32(42);
  t('rng: seed が同じなら同じ列', r1() === r2() && r1() === r2());
  t('rng: seed が違えば違う列', mulberry32(1)() !== mulberry32(2)());

  // --- permP
  const pa = permP([1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], 5000, 7);
  t('permP: 完全逆順は最小 p(2/(n!)=1/360 付近)', pa > 0 && pa < 0.01);
  t('permP: 決定的(2 回同じ)',
    permP([1, 2, 3, 4, 5], [2, 1, 4, 3, 5], 3000, 11) === permP([1, 2, 3, 4, 5], [2, 1, 4, 3, 5], 3000, 11));
  t('permP: p は 0 にならない', permP([1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7], 2000, 3) > 0);
  {
    // 無相関に近い列では p が大きいこと(p ≧ 0.2)
    const pb = permP([1, 2, 3, 4, 5, 6, 7, 8], [3, 1, 4, 8, 2, 7, 5, 6], 5000, 5);
    t('permP: ばらばらな列は p が大きい', pb > 0.2);
  }

  // --- holm
  const h = holm([0.01, 0.04, 0.03]);
  // 昇順 0.01(×3=0.03), 0.03(×2=0.06), 0.04(×1=0.04→単調化で 0.06)
  t('holm: 最小 p の調整', close(h[0], 0.03, 1e-12));
  t('holm: 単調化(0.04 が 0.06 に持ち上がる)', close(h[1], 0.06, 1e-12));
  t('holm: 中間', close(h[2], 0.06, 1e-12));
  t('holm: 1 で頭打ち', holm([0.9, 0.9, 0.9]).every(v => v === 1));
  t('holm: m=1 は素通し', close(holm([0.02])[0], 0.02, 1e-12));

  // --- partial
  t('partial: 交絡ゼロなら素通し', close(partial(0.5, 0, 0), 0.5, 1e-12));
  t('partial: 完全交絡で消える', close(partial(0.5, 0.5, 1.0), 0, 1e-9) || Number.isNaN(partial(0.5, 0.5, 1.0)));

  // --- icc1
  const same = icc1([[1, 1], [5, 5], [9, 9]]);   // 束内分散 0 ⇒ ICC = 1
  t('icc: 束内が同一なら ICC = 1', close(same.icc, 1, 1e-9));
  const none = icc1([[1, 9], [1, 9], [1, 9]]);   // 束間分散 0 ⇒ ICC ≦ 0
  t('icc: 束間に差が無ければ ICC ≦ 0', none.icc <= 1e-9);
  t('icc: DEFF は 1 以上', same.deff >= 1 && none.deff >= 1);
  t('icc: 束 1 つでは NaN', Number.isNaN(icc1([[1, 2, 3]]).icc));

  // --- bundleKey
  t('bundle: 枝番を落とす(段1c → 段)', bundleKey('段1c 形式群則の混合項と hρ') === '段');
  t('bundle: Λ9 → Λ', bundleKey('Λ9 円分子を副有限捩れから') === 'Λ');
  t('bundle: B1 → B', bundleKey('B1 アーベル閉包の定義') === 'B');
  t('bundle: 数字が無ければそのまま', bundleKey('hρ を外す帳簿 3 手') === 'hρ');
  t('bundle: B2+B3 は先頭側で束ねる', bundleKey('B2+B3 E_σ と K^ab の分解') === 'B');

  // --- ★COST(自己申告)の書式
  t('cost: 素の 1 行', JSON.stringify(parseCostLine('COST[Y19e]: 安')) === '{"key":"Y19e","level":"安","slot":null,"note":""}');
  t('cost: ★や - の飾りを許す', parseCostLine('- ★COST[B2+B3]: 高 — 抽象核が 2 本増えた')?.level === '高');
  t('cost: 一言を切り出す', parseCostLine('COST[Λ12]: 並 — 見積どおり')?.note === '見積どおり');
  t('cost: 全角コロン', parseCostLine('COST[Y9]：安')?.key === 'Y9');
  t('cost: 書式の見本は数えない(値も見本)', parseCostLine('COST[<持ち場>]: <安|並|高>') === null);
  // ★★値が正しくて**鍵だけが見本**の行。★これが無いと `[<>]` の番人を外しても鳴らない
  //   (第 24 回がわざと壊して気づいた。値の側の列挙で弾かれていて空虚だった)。
  t('cost: ★鍵だけが見本でも数えない', parseCostLine('COST[<持ち場>]: 安 — 見本') === null);
  t('cost: 値が列挙外なら読まない', parseCostLine('COST[X]: そこそこ') === null);
  t('cost: COST を含むだけの散文は読まない', parseCostLine('この行は COST[ を含むが書式ではない') === null);
  t('cost: 鍵の前後の空白を落とす', parseCostLine('COST[ Y19e ]: 高')?.key === 'Y19e');
  t('cost: 橋が無ければ slot は null', parseCostLine('COST[Y19e]: 安 — x')?.slot === null);
  t('cost: ★明示の橋を読む', parseCostLine('COST[LT]: 安 | 持ち場=局所Tate双対性の道')?.slot === '局所Tate双対性の道');
  t('cost: 橋を一言から取り除く',
    parseCostLine('COST[LT]: 安 | 持ち場=局所Tate双対性の道 — 早かった')?.note === '早かった');
  {
    const tmp3 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest3.md';
    fs.writeFileSync(tmp3, ['COST[LT]: 安 | 持ち場=局所Tate双対性の道'].join('\n'));
    const rows4 = [{ description: '局所Tate双対性の道', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } }];
    const buf4 = []; const real4 = console.log;
    console.log = (...a) => buf4.push(a.join(' '));
    try { cmdCost(rows4, [tmp3]); } finally { console.log = real4; }
    t('cost: ★橋があれば鍵が違っても当たる', buf4.join('\n').includes('突き合わせ 1 件'));
    fs.unlinkSync(tmp3);
  }
  // ★突き合わせの境界(ここが緩いと Y19 が Y19e を食う)
  t('join: 完全一致', matchesKey('Y19e', 'Y19e'));
  t('join: 鍵 + 空白 + 説明', matchesKey('Y19e compat の心臓', 'Y19e'));
  t('join: ★Y19 は Y19e を拾わない', !matchesKey('Y19e compat の心臓', 'Y19'));
  t('join: Y19 は Y19 の持ち場を拾う', matchesKey('Y19 絶対Galois群の分岐', 'Y19'));
  t('join: 記号を含む鍵', matchesKey('B2+B3 E_σ と K^ab の分解', 'B2+B3'));
  t('join: 空の鍵は拾わない', !matchesKey('Y19e', ''));
  {
    const tmp = process.env.TEMP || '.';
    const p = tmp + '/_agent-timing-cost-selftest.md';
    fs.writeFileSync(p, ['# 見本', 'COST[<持ち場>]: <安|並|高>', 'COST[Y19e]: 安 — 早かった',
                         '★COST[Λ12]: 高', 'COST[Y19e]: 並 — 書き直した'].join('\n'));
    const { costs } = readCosts([p]);
    t('readCosts: 見本を除いて 2 鍵', costs.size === 2);
    t('readCosts: 同じ鍵は後の記述で上書き', costs.get('Y19e').level === '並');
    t('readCosts: 行番号を持つ', costs.get('Λ12').line === 4);
    const { missing } = readCosts([p + '.nope']);
    t('readCosts: 無い記録は missing に', missing.length === 1);
    fs.unlinkSync(p);
  }
  {
    // ★★空虚さの見張り: 申告が 0 件のとき **判定めいた文を出さない**
    const grab2 = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const out0 = grab2(() => cmdCost([{ description: 'Y19e x', durationMs: 1, lines: 1, est: null }], ['/nope/none.md']));
    t('cost: 申告 0 件なら「何も言わない」と書く', out0.includes('分母が立たない'));
    t('cost: 申告 0 件で ρ を出さない', !out0.includes('Spearman'));
    const tmp2 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest2.md';
    fs.writeFileSync(tmp2, ['COST[Y1]: 安', 'COST[Y2]: 高', 'COST[ZZ]: 並'].join('\n'));
    const rows0 = [{ description: 'Y1 なにか', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } },
                   { description: 'Y2 べつ', durationMs: 60000, lines: 300, est: { lo: 80, hi: 120 } }];
    const out1 = grab2(() => cmdCost(rows0, [tmp2]));
    t('cost: 突き合わせ件数を出す', out1.includes('突き合わせ 2 件'));
    t('cost: 当たらなかった申告を名指しする', out1.includes('ZZ'));
    t('cost: ★件数不足なら判定を出さない', out1.includes('判定を出さない') && !out1.includes('Spearman'));
    // ★★水準が偏っていたら、件数が足りていても判定を出さない(事前登録した (b))
    {
      const tmp4 = (process.env.TEMP || '.') + '/_agent-timing-cost-selftest4.md';
      const keys = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'B1'];
      fs.writeFileSync(tmp4, keys.map((k, i) => `COST[${k}]: ${i < 7 ? '安' : '並'}`).join('\n'));
      const rows4 = keys.map(k => ({ description: k + ' x', durationMs: 60000, lines: 100, est: { lo: 80, hi: 120 } }));
      const out4 = grab2(() => cmdCost(rows4, [tmp4]));
      t('cost: ★8 件あっても水準が偏れば判定を出さない', out4.includes('突き合わせ 8 件')
        && out4.includes('判定を出さない') && !out4.includes('Spearman'));
    }
    fs.unlinkSync(tmp2);
  }

  // --- ★族(事前登録)—— 欄が黙って消えたら鳴る
  t('family: 事前登録した欄が族に居る(estMid)', FAMILY.some(f => f.key === 'estMid'));
  t('family: 族は 7 本', FAMILY.length === 7);
  t('family: 鍵が重複しない', new Set(FAMILY.map(f => f.key)).size === FAMILY.length);
  // ★★メタ第 25 回が決めた正典を固定する。黙って差し替わったら鳴る(理由は FAMILY の頭に書いた)。
  t('family: ★行数の正典は lines(linesWritten ではない)', FAMILY.some(f => f.key === 'lines'));
  t('family: ★linesWritten は族に入れない(m を増やさない)', !FAMILY.some(f => f.key === 'linesWritten'));
  // ★★COST(自己申告)は n が足りるまで族に入れない —— 事前登録(メタ第 25 回)。
  //    条件は 2 つとも満たすこと: (a) 突き合わせ 8 件以上 (b) 各水準に 3 件以上。
  t('family: ★COST はまだ族に入れない', !FAMILY.some(f => f.key === 'costLevel'));

  // ── M142(メタ第 30 回): heredoc を剥ぐ / 分母の突然変異 ────────────────
  t('heredoc: 素の本文は残る', stripHeredocs('lake build ABC3').includes('lake build'));
  t('heredoc: ★本文の中の lake build が消える',
    !stripHeredocs("cat > a.md <<'EOF'\nlake build ABC3\nEOF\ngit status").includes('lake build'));
  t('heredoc: ★heredoc の外の git は残る',
    stripHeredocs("cat > a.md <<'EOF'\nlake build ABC3\nEOF\ngit status").includes('git status'));
  t('heredoc: 引用符なしの tag も剥ぐ',
    !stripHeredocs('cat > a.md <<EOF\nlake build\nEOF').includes('lake build'));
  // ★引用符なしの tag を取り落とすと「終端が来ない」ことになり、★後ろが丸ごと消える(黙って分母が減る)
  t('heredoc: ★引用符なしでも終端の後ろが残る',
    stripHeredocs('cat > a.md <<EOF\nlake build\nEOF\ngit status').includes('git status'));
  t('heredoc: <<- も剥ぐ', !stripHeredocs('cat <<-EOF\nlake build\nEOF').includes('lake build'));
  t('heredoc: 閉じないまま終わっても落ちない', stripHeredocs("cat <<'EOF'\nlake build") === 'cat <<\'EOF\'');
  t('heredoc: 空/未定義でも落ちない', stripHeredocs(undefined) === '' && stripHeredocs('') === '');
  t('heredoc: ★分類が変わる(素は lake build、剥ぐと git)',
    classifyBash("git commit -F- <<'EOF'\nlake build ABC3\nEOF") === 'lake build'
    && classifyBash(stripHeredocs("git commit -F- <<'EOF'\nlake build ABC3\nEOF")) === 'git');
  t('heredoc: ★既定の classifyBash は変えていない', classifyBash("x <<'E'\nlake build\nE") === 'lake build');
  {
    // trapMutate: Bash の命令文だけを置き換え、他の行はそのまま返す
    const l1 = JSON.stringify({ message: { content: [{ type: 'tool_use', name: 'Bash', id: 'a', input: { command: 'ls' } }] } });
    const l2 = JSON.stringify({ message: { content: [{ type: 'tool_use', name: 'Read', id: 'b', input: { file_path: 'x.lean' } }] } });
    const m = trapMutate([l1, l2, 'not json'].join('\n'), 'ZZZ');
    t('trap: Bash の件数を数える', m.n === 1);
    t('trap: ★命令文が置き換わる', m.text.includes('ZZZ'));
    t('trap: Bash 以外の行は素通り', m.text.includes('x.lean') && m.text.includes('not json'));
  }
  {
    // covariatesFromText の opts が既定では現行と同じで、明示すると効くこと
    const mk = (name) => JSON.stringify({ message: { content: [{ type: 'tool_use', name, id: name, input: { file_path: 'D:/r/lean/ABC3/Found/A.lean' } }] } });
    const chk = [
      JSON.stringify({ message: { content: [{ type: 'tool_use', name: 'mcp__x__lean_check', id: 'c1', input: {} }] } }),
      JSON.stringify({ message: { content: [{ type: 'tool_result', tool_use_id: 'c1', content: 'OK (0.1 秒)\ninfo: Foo.error_bound' }] } }),
    ];
    const text = [mk('Read'), mk('Read'), mk('Write'), ...chk].join('\n');
    const a = covariatesFromText(text, '.');
    const b = covariatesFromText(text, '.', { fileSource: 'written' });
    t('cov: 既定は Read も数える(files に載る)', a.files.length === 1 && a.leanChecks === 1);
    t('cov: written に限っても Write があれば載る', b.files.length === 1);
    t('cov: ★字面規則では OK でも error を含めば失敗に数える', a.checkFails === 1);
    t('cov: ★見出し規則なら失敗に数えない', covariatesFromText(text, '.', { failRule: 'header' }).checkFails === 0);
    const textRO = [mk('Read'), mk('Read')].join('\n');
    t('cov: ★読んだだけの本は written では落ちる',
      covariatesFromText(textRO, '.').files.length === 1
      && covariatesFromText(textRO, '.', { fileSource: 'written' }).files.length === 0);
    const textME = [mk('Read'), mk('MultiEdit')].join('\n');
    t('cov: ★MultiEdit も「自分が書いた」に数える',
      covariatesFromText(textME, '.', { fileSource: 'written' }).files.length === 1);
  }
  {
    // familyLadder は cmdExplain と同じ規則(MIN_N 以上・全件同値でない)
    const rows = Array.from({ length: 20 }, (_, i) => ({
      durationMs: (i + 1) * 1000, lines: (i + 1) * 10, toolUses: i + 1, tokens: (i + 1) * 7,
      cores: i % 3, leanChecks: i % 4, checkFails: i % 5, estMid: 100,
    }));
    const L = familyLadder(rows);
    t('ladder: 全件同値の欄は落ちる(estMid)', !L.some(x => x.key === 'estMid'));
    t('ladder: 6 本が残る', L.length === 6);
    t('ladder: Holm は p 以上', L.every(x => x.holm >= x.p - 1e-12));
    t('ladder: ★Holm は少なくとも 1 本で p を厳密に上回る(補正が効いている)',
      L.some(x => x.holm > x.p + 1e-9));
    t('ladder: 単調(lines)は言える', L.find(x => x.key === 'lines').say === true);
    t('ladder: ★言えない欄がある(閾が効いている)', L.some(x => x.say === false));
    t('ladder: n が足りない欄は落ちる', familyLadder(rows.slice(0, 5)).length === 0);
  }

  // ── ★★★M149 の事前登録(メタ第 31 回)。★黙って書き換わったら鳴る ───────────
  {
    t('M149: 締切は 2026-09-07T16:35:29Z のまま', M149_CUTOFF === '2026-09-07T16:35:29Z');
    t('M149: 必要件数は 134 のまま', M149_NEED === 134);
    t('M149: 一次の欄は cores', M149_PRIMARY === 'cores');
    t('M149: 分母は written(touched ではない)', M149_DENOM === 'written');
    // ★導出を式で固定する(定数だけ直しても鳴るように)
    const need = (rho, alpha) => {
      // Acklam の逆正規(±4.5e-4 の精度。★閾の 134 を判定するには十分)
      const inv = (p) => {
        const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
        const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01];
        const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
          -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
        const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
        const pl = 0.02425;
        if (p < pl) { const q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
        if (p > 1 - pl) return -inv(1 - p);
        const q = p - 0.5, r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
          (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
      };
      return 3 + ((inv(1 - alpha / 2) + inv(0.80)) / Math.atanh(rho)) ** 2;
    };
    t('M149: ★必要件数 134 は「ρ=0.30・検出力 0.80・α=0.05/7」から出る',
      Math.ceil(need(0.30, 0.05 / FAMILY.length)) === M149_NEED);
    t('M149: ★ρ が小さいほど必要件数は増える', need(0.20, 0.05 / 7) > need(0.30, 0.05 / 7));
    t('M149: ★α が緩いほど必要件数は減る', need(0.30, 0.05) < need(0.30, 0.05 / 7));
    // freshRows —— 締切の前後で確実に割れる
    const rr = [{ ts: '2026-09-07T16:35:28Z' }, { ts: '2026-09-07T16:35:29Z' },
      { ts: '2026-09-07T16:35:30Z' }, { ts: null }, {}];
    t('M149: freshRows は締切より後だけ', freshRows(rr).length === 1
      && freshRows(rr)[0].ts === '2026-09-07T16:35:30Z');
    t('M149: freshRows は ts が無い行を落とす', freshRows([{ ts: null }, {}]).length === 0);
    // ★filterByTime(第 31 回の突然変異 A10 が素通りしたので純関数に割って試験を足した)
    t('M149: --since は境界を含まない', filterByTime(rr, '2026-09-07T16:35:29Z').length === 1);
    t('M149: --until は境界を含む', filterByTime(rr, null, '2026-09-07T16:35:29Z').length === 2);
    t('M149: --since と --until を両方効かせる',
      filterByTime(rr, '2026-09-07T16:35:28Z', '2026-09-07T16:35:29Z').length === 1);
    t('M149: どちらも無ければ素通し', filterByTime(rr).length === rr.length);
    // ★停止規則 —— 届かなければ ladder を **作らない**
    const mk = (n) => Array.from({ length: n }, (_, i) => ({
      durationMs: (i + 1) * 1000, lines: i + 1, cores: i % 7, toolUses: i + 2,
      tokens: (i + 1) * 10, leanChecks: i % 5, checkFails: i % 3, estMid: 100 + i,
    }));
    const gLo = m149Gate(mk(M149_NEED - 1));
    t('M149: ★届かなければ open=false', gLo.open === false);
    t('M149: ★届かなければ ladder は null(計算すらしない)', gLo.ladder === null);
    t('M149: ★不足数を出す', gLo.short === 1 && gLo.n === M149_NEED - 1);
    const gHi = m149Gate(mk(M149_NEED));
    t('M149: ★届けば open=true', gHi.open === true);
    t('M149: ★届けば ladder は族の階段', Array.isArray(gHi.ladder) && gHi.ladder.length >= 1
      && gHi.ladder.every(x => Number.isFinite(x.holm)));
    t('M149: ★届いたときの不足数は 0', gHi.short === 0);
    // ★分母が written であることを実データ形式で確かめる(touched なら Read だけの本が代表になる)
    {
      const mkLine = (o) => JSON.stringify({ message: { role: 'assistant', content: o } });
      const rd = { type: 'tool_use', id: 'r1', name: 'Read', input: { file_path: 'lean/ABC3/Found/A.lean' } };
      const wr = { type: 'tool_use', id: 'w1', name: 'Write', input: { file_path: 'lean/ABC3/Found/B.lean', content: 'x\ny\n' } };
      const txt = [mkLine([rd]), mkLine([rd]), mkLine([rd]), mkLine([wr])].join('\n');
      const cT = covariatesFromText(txt, '.', {});
      const cW = covariatesFromText(txt, '.', { fileSource: 'written' });
      t('M149: touched は Read だけの本を代表にする(これが M145 の誤り)', cT.file === 'lean/ABC3/Found/A.lean');
      t('M149: ★written は自分が書いた本を代表にする', cW.file === 'lean/ABC3/Found/B.lean');
      t('M149: ★written は Read しかしていない agent で欠測する',
        covariatesFromText([mkLine([rd])].join('\n'), '.', { fileSource: 'written' }).file === null);
    }
    // ★M111 は族に入れない(M149 手順 3 の決着を固定する)
    t('M149: ★M111(完了時点の行数)は族に入れない', !FAMILY.some(f => /^lines(Written|AtDone)$/.test(f.key)));
  }
  t('cost: ★族に入れる条件(8 件かつ各水準 3 件)', costFamilyReady([]) === false
     && costFamilyReady(['安', '安', '安', '並', '並', '並', '高', '高']) === false
     && costFamilyReady(['安', '安', '安', '並', '並', '並', '高', '高', '高']) === true);

  // --- ★欄ごとの complete-case
  {
    const rs = [{ a: 1, b: NaN }, { a: 2, b: 5 }, { a: 3, b: NaN }];
    t('usableRows: NaN の行を落とす', usableRows(rs, 'b').length === 1 && usableRows(rs, 'a').length === 3);
    t('MIN_N は 8', MIN_N === 8);
  }
  {
    // ★★空虚さの見張り: 欠測のある欄を足しても **他の欄の ρ が動かない** こと。
    //   ★ここが壊れる(listwise に戻す/欠測行を 0 で埋める)と、この 2 件が鳴る。
    const mk = (n) => Array.from({ length: n }, (_, i) => ({
      durationMs: (i + 1) * 1000, lines: (i + 1) * 10, toolUses: i + 1, tokens: (i + 1) * 7,
      cores: i % 3, leanChecks: i % 4, checkFails: i % 5,
      // ★半分だけ見積がある(欠測が欄で違う状況を作る)
      estMid: i % 2 === 0 ? 100 * (n - i) : NaN,
      description: 'X' + i,
    }));
    const grab = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const rows = mk(20);
    const out = grab(() => cmdExplain(rows));
    t('explain: 欠測のある欄も検定に入る(m = 7)', out.includes('m = 7'));
    t('explain: 欄ごとの n を印字する(10 と 20 が並ぶ)', /\s10\s/.test(out) && /\s20\s/.test(out));
    // 同じ盤面から estMid を丸ごと外したものと、lines の ρ を比べる
    const rows2 = mk(20).map(r => ({ ...r, estMid: NaN }));
    const rho1 = spearman(rows.map(r => r.lines), rows.map(r => r.durationMs));
    const rho2 = spearman(rows2.map(r => r.lines), rows2.map(r => r.durationMs));
    t('explain: 欠測欄を足しても他の欄の ρ は動かない', close(rho1, rho2, 1e-12));
    const out2 = grab(() => cmdExplain(rows2));
    t('explain: 全欠測の欄は「測れなかった」に落ちる', out2.includes('測れなかった') && out2.includes('m = 6'));
    // ★件数不足(MIN_N 未満)の欄はそう書く
    const rows3 = mk(20).map((r, i) => ({ ...r, estMid: i < 5 ? r.estMid : NaN }));
    const out3 = grab(() => cmdExplain(rows3));
    t('explain: MIN_N 未満は件数不足と書いて m から外す', out3.includes('件数不足') && out3.includes('m = 6'));
    t('explain: 向きを断定する文を出さない', !/(長いほど|多いほど|ほど遅い|ほど速い)/.test(out));
  }

  // --- describe
  const d = describe([1, 2, 3, 4]);
  t('describe: 中央値 2.5', close(d.med, 2.5));
  t('describe: 合計 10', d.sum === 10);

  // --- 5b. --diag(診断の族 × 費用。★メタ第 27 回)
  {
    // Fisher —— 既知の値で較正する
    t('diag: Fisher お茶の実験 [[3,1],[1,3]] = 0.4857', close(fisher(3, 1, 1, 3), 0.4857142857, 1e-6));
    t('diag: Fisher 完全分離 [[10,0],[0,10]]', close(fisher(10, 0, 0, 10), 2 / 184756, 1e-9));
    t('diag: Fisher 独立な表は 1', close(fisher(5, 5, 5, 5), 1, 1e-9));
    t('diag: Fisher 空の行は 1', fisher(0, 0, 3, 3) === 1);
    t('diag: Fisher は対称', close(fisher(7, 2, 3, 8), fisher(2, 7, 8, 3), 1e-12));
    // needBalancedN —— 割合が同じなら null、違えば α を割る n を返す
    t('diag: あと何件 割合が同じなら null', needBalancedN(5, 5, 5, 5, 0.005) === null);
    const nb = needBalancedN(8, 2, 2, 8, 0.005);
    t('diag: あと何件 は MIN_N 以上の整数', Number.isInteger(nb) && nb >= MIN_N);
    t('diag: あと何件 は本当に α を割る', (() => { const A = Math.round(0.8 * nb), C = Math.round(0.2 * nb); return fisher(A, nb - A, C, nb - C) <= 0.005; })());
    t('diag: あと何件 その 1 つ手前では割らない', (() => { const n = nb - 1; if (n < MIN_N) return true; const A = Math.round(0.8 * n), C = Math.round(0.2 * n); return fisher(A, n - A, C, n - C) > 0.005; })());
    // tallyDiagnostics —— doc / meta / 本体セッション(ag 空)は入れない
    const errs = [
      { ag: 'T1', src: 'tree', msg: 'failed to synthesize instance Foo' },
      { ag: 'T1', src: 'scratch', msg: 'unsolved goals x' },
      { ag: 'T1', src: 'doc', msg: 'failed to synthesize (文書を印字しただけ)' },
      { ag: 'T1', src: 'meta', msg: 'failed to synthesize (改善係の probe)' },
      { ag: '', src: 'tree', msg: 'failed to synthesize (本体セッション)' },
      { ag: 'T2', src: 'scratch', msg: 'Type mismatch here' },
    ];
    const tal = tallyDiagnostics(errs);
    t('diag: 本体セッション(ag 空)は数えない', !tal.has(''));
    t('diag: doc / meta は数えない', tal.get('T1').B.total === 2);
    // ★スコープが挙げても doc / meta は数えない(★門番が二重であることを確かめる)
    const talX = tallyDiagnostics(errs, DIAG_FAMILIES, [{ id: 'X', srcs: ['tree', 'scratch', 'doc', 'meta'] }]);
    t('diag: スコープが挙げても doc / meta は数えない', talX.get('T1').X.total === 2);
    t('diag: スコープ A は tree だけ', tal.get('T1').A.total === 1 && tal.get('T1').A.fam[0] === 1 && tal.get('T1').A.fam[1] === 0);
    t('diag: スコープ B は tree + scratch', tal.get('T1').B.fam[0] === 1 && tal.get('T1').B.fam[1] === 1);
    t('diag: 族は部分文字列で当てる', tal.get('T2').B.fam[4] === 1);
    // medianSplit —— 同値は low 側(★事前に決めた向き)
    const ms2 = medianSplit([1, 1, 1, 5]);
    t('diag: 中央値 2 値化 同値は low', ms2.isHigh(1) === false && ms2.isHigh(5) === true);
    // cmdDiag —— 件数不足なら検定しない / 向きを断定しない / m を印字する
    const gr = (fn) => { const buf = []; const real = console.log;
      console.log = (...a) => buf.push(a.join(' ')); try { fn(); } finally { console.log = real; }
      return buf.join('\n'); };
    const mkRec = (i, fam) => ({ toolUseId: 'A' + i, agentType: DIAG_TYPE, toolUses: 10, durationMs: (i + 1) * 60000, day: '2026-09-0' + (1 + (i % 3)), fam });
    const recsA = Array.from({ length: 24 }, (_, i) => mkRec(i, i % 2 === 0));
    const digA = { v: 3, errors: recsA.filter(r => r.fam).map(r => ({ ag: r.toolUseId, src: 'tree', msg: 'failed to synthesize instance' })) };
    const outD = gr(() => cmdDiag(recsA, digA));
    t('diag: m = 10 を印字する', outD.includes('m = 10'));
    t('diag: 曝露が MIN_N を満たせば検定する', /言える|言えない/.test(outD));
    t('diag: 曝露ゼロの族は件数不足と書く', outD.includes('件数不足'));
    t('diag: 向きを断定する文を出さない', !/(ほど遅い|ほど速い|のほうが高い|のほうが遅い|を食っている)/.test(outD));
    t('diag: 「族が時間を食う」は否定形でしか書かない', /が時間を食う」とは言えない/.test(outD));
    // ★★多重比較の補正が本当に効いているか —— 素の p は 0.05 を割るが、m = 10 を掛けると割らない盤面。
    //   ★ここが素通りすると「言えない」が「★言える」に化ける。
    const hi = new Set([...Array.from({ length: 9 }, (_, k) => 12 + k), 0, 1, 2]);   // 曝露 12 件(うち高 9)
    const recsC = Array.from({ length: 24 }, (_, i) => mkRec(i, hi.has(i)));
    const digC = { v: 3, errors: recsC.filter(r => r.fam).map(r => ({ ag: r.toolUseId, src: 'tree', msg: 'failed to synthesize instance' })) };
    const outC = gr(() => cmdDiag(recsC, digC));
    t('diag: 盤面が意図どおり(素の p ≒ 0.039)', /0\.039/.test(outC));
    // ★「★言えるのは…」という地の文にも同じ字面が入るので、★行末の判定欄だけを見る。
    t('diag: 補正すると言えない(補正を外すと言えるになる盤面)', !/★言える\s*$/m.test(outC));
    const digEmpty = { v: 3, errors: [] };
    const outE = gr(() => cmdDiag(recsA, digEmpty));
    t('diag: 診断ゼロでも落ちない(全部 件数不足)', (outE.match(/件数不足/g) || []).length === 10);
    t('diag: 突き合わせられない件数を印字する', outD.includes('本体セッション'));
  }

  // --- 5c. --main(本体セッションの費用。★事前登録した単位が本当にその単位か)
  {
    t('main: normLean 絶対パス', normLean('D:\\Math_ABC3\\lean\\ABC3\\Found\\X.lean') === 'ABC3/Found/X.lean');
    t('main: normLean 相対パス', normLean('lean/ABC3/Found/PGC/ArtinMap.lean') === 'ABC3/Found/PGC/ArtinMap.lean');
    t('main: normLean ABC3 が無ければそのまま', normLean('tools/check.mjs') === 'tools/check.mjs');
    t('main: leanTarget は file_path から', leanTarget('Edit', { file_path: 'D:/Math_ABC3/lean/ABC3/Found/A.lean' }) === 'ABC3/Found/A.lean');
    t('main: leanTarget は .lean 以外の file_path を採らない', leanTarget('Read', { file_path: 'ResearchPaper/x.md' }) === null);
    t('main: leanTarget は Bash の本文からも拾う', leanTarget('Bash', { command: 'lake build ABC3.Found.B 2>&1 | head; wc -l lean/ABC3/Found/B.lean' }) === 'ABC3/Found/B.lean');
    t('main: leanTarget 無ければ null', leanTarget('Bash', { command: 'git status' }) === null);
    t('main: classifyBash は lake build を node より先に当てる', classifyBash('cd x && lake build ABC3 && node tools/check.mjs') === 'lake build');
    t('main: classifyBash 既定は その他', classifyBash('ls -la') === 'その他');
    // bundleByFile —— ファイルで分かれ、gap で切れる
    const mk = (id, t0, ms, tgt, n = 'Bash', bk = '') => ({ id, s: 'S', n, bk, t0, t1: t0 + ms, tgt, e: 0 });
    const cs = [mk('1', 0, 1000, 'ABC3/A.lean'), mk('2', 60000, 1000, 'ABC3/A.lean'),
                mk('3', 120000, 1000, 'ABC3/B.lean'),
                mk('4', 60000 * 100, 1000, 'ABC3/A.lean'), mk('5', 0, 1, null)];
    const bs = bundleByFile(cs, 30);
    t('main: 塊はファイルで分かれる', bs.filter(b => b.file === 'ABC3/B.lean').length === 1);
    t('main: 30 分空いたら別の塊', bs.filter(b => b.file === 'ABC3/A.lean').length === 2);
    t('main: 対象の無い呼び出しは塊に入らない', bs.reduce((s, b) => s + b.calls.length, 0) === 4);
    // delegLike —— (a)(b)(c) が独立に効く
    const many = (file, k, lean) => ({ file, start: 0, end: 1000 * k, calls: Array.from({ length: k }, (_, i) =>
      mk('x' + i, i * 1000, 100, file, lean && i === 0 ? 'mcp__abc3-lean__lean_check' : 'Bash')) });
    t('main: (a) 短い塊は輪郭に入らない', delegLike(many('ABC3/A.lean', 5, true), [], 20).ok === false);
    t('main: (c) Lean を叩かない塊は入らない', delegLike(many('ABC3/A.lean', 30, false), [], 20).why === 'Lean を叩いていない');
    t('main: 3 条件を満たせば入る', delegLike(many('ABC3/A.lean', 30, true), [], 20).ok === true);
    const foreign = [{ ...mk('w', 5000, 10, 'ABC3/OTHER.lean', 'Edit'), }];
    t('main: (b) 期間内に他の .lean を書いたら入らない', /他の \.lean/.test(delegLike(many('ABC3/A.lean', 30, true), foreign, 20).why));
    t('main: (b) 同じファイルへの Write は落とさない',
      delegLike(many('ABC3/A.lean', 30, true), [{ ...mk('w', 5000, 10, 'ABC3/A.lean', 'Edit') }], 20).ok === true);
    // activeMs —— 頭打ちが効く
    // ★t0 = 0 は「時刻が取れなかった」印なので activeMs は外す。★試験も 0 を使わない。
    t('main: 稼働の代理は間隔を cap で頭打ちにする',
      activeMs([mk('a', 1000, 1000, null), mk('b', 3601000, 1000, null)], 300) === 1000 + 300000 + 1000);
    t('main: 時刻の無い呼び出し(t0=0)は稼働に数えない', activeMs([mk('a', 0, 1000, null)], 300) === 0);
    // attributeDiagnostics —— 秒で割る / 割れないものは黙って寄せない / agent の分は扱わない
    const cA = mk('c1', 10000, 2000, 'ABC3/A.lean');            // t1 = 12000 → 秒 12
    const cB = mk('c2', 11000, 1000, 'ABC3/B.lean');            // t1 = 12000 → 秒 12
    const E = (ts, file, msg, ag = '') => ({ ts, file, msg, src: 'tree', ag });
    const a1 = attributeDiagnostics([cA], [E(12, 'lean/ABC3/A.lean', 'unsolved goals')]);
    t('main: 診断を秒で呼び出しに割り当てる', a1.matched === 1 && a1.byCall.get('c1').fam[1] === 1);
    const a2 = attributeDiagnostics([cA, cB], [E(12, 'lean/ABC3/B.lean', 'unsolved goals')]);
    t('main: 同じ秒でもファイルで割れる', a2.matched === 1 && a2.byCall.has('c2'));
    const a3 = attributeDiagnostics([cA, cB], [E(12, 'lean/ABC3/Z.lean', 'unsolved goals')]);
    t('main: 割れないものは ambiguous に数える(寄せない)', a3.matched === 0 && a3.ambiguous === 1);
    const a4 = attributeDiagnostics([cA], [E(12, 'lean/ABC3/A.lean', 'unsolved goals', 'toolu_agent')]);
    t('main: agent の診断は本体に数えない', a4.matched === 0);
    const a5 = attributeDiagnostics([cA], [E(99, 'lean/ABC3/A.lean', 'unsolved goals')]);
    t('main: 対応する呼び出しが無ければ unmatched', a5.unmatched === 1);
    const a6 = attributeDiagnostics([cA], [{ ts: 12, file: 'x', msg: 'unsolved goals', src: 'scratch', ag: '' }]);
    t('main: 既定のスコープは tree だけ', a6.matched === 0);
  }

  // ────────────────────────────────────────────────────────────
  // ★★M159 —— 「他と混ざらない命令だけ」の規則(SOLO-1)の較正
  // ★★第 30 回・第 31 回の反省を踏まえ、**素通りしそうな側を先に書く**:
  //   混ざっている命令を落とせるか / heredoc の幽霊を落とせるか / 時間の付け方。
  // ────────────────────────────────────────────────────────────
  {
    const NL = String.fromCharCode(10);
    t('solo: 単独の check.mjs は当たる', soloOf('node tools/check.mjs --brief') === 'check.mjs');
    t('solo: Windows の \\ 区切りでも当たる', soloOf('node tools\\check.mjs --brief') === 'check.mjs');
    t('solo: ★2 本混ざったら落とす',
      soloOf('node tools/check.mjs --brief && node tools/graph.mjs') === null);
    t('solo: ★lake build が混ざったら落とす',
      soloOf('lake build ABC3 && node tools/check.mjs') === null);
    t('solo: ★lake env lean が混ざったら落とす',
      soloOf('lake env lean X.lean; node tools/graph.mjs') === null);
    t('solo: どれにも当たらなければ null', soloOf('git status') === null);
    t('solo: 空でも落ちない', soloOf('') === null && soloOf(null) === null && soloOf(undefined) === null);
    t('solo: ★git や grep が繋がっていても L1 では単独(★既知の穴。文書化済)',
      soloOf('git status && node tools/check.mjs --brief') === 'check.mjs');
    // ★heredoc の幽霊(M146)。本文に別の道具名が書き写されている記録。
    {
      const cmd = ['cat > x.md <<EOF', 'node tools/graph.mjs をここで叩いた', 'EOF',
        'node tools/check.mjs --brief'].join(NL);
      t('solo: ★L0(素)は heredoc の本文に騙されて落とす', soloOf(cmd, 0) === null);
      t('solo: ★★L1 は heredoc を剥いで正しく check.mjs に当てる', soloOf(cmd, 1) === 'check.mjs');
    }
    // ★L2(厳格)—— 他の tools/ を許さない
    t('solo: L1 は他の tools/ を許す',
      soloOf('node tools/check.mjs && node tools/hedge-index.mjs', 1) === 'check.mjs');
    t('solo: ★L2 は他の tools/ を許さない',
      soloOf('node tools/check.mjs && node tools/hedge-index.mjs', 2) === null);
    t('solo: L2 でも本当に単独なら当たる', soloOf('node tools/check.mjs --brief', 2) === 'check.mjs');
    t('solo: ★L2 は python が混ざったら落とす',
      soloOf('node tools/check.mjs && python x.py', 2) === null
      && soloOf('node tools/check.mjs && python x.py', 1) === 'check.mjs');
    t('solo: 同じ道具を 2 回書いても 1 件(★既知の穴。文書化済)',
      soloOf('node tools/check.mjs; node tools/check.mjs') === 'check.mjs');
    t('solo: mentioned は基底名を拾う',
      [...soloToolsMentioned('node tools/a.mjs && node tools\\b.py')].sort().join(',') === 'a.mjs,b.py');
    // ★集計(tallySolo)—— 時間の付け方と、母集団から落ちるもの
    {
      const C = (cmd, t0, t1) => ({ cmd, t0, t1 });
      const cs = [
        C('node tools/check.mjs', 1000, 8000),          // 7 秒
        C('node tools/check.mjs', 1000, 4000),          // 3 秒
        C('node tools/graph.mjs', 1000, 4000),          // 3 秒
        C('node tools/check.mjs && node tools/graph.mjs', 1000, 100000),  // ★混ざり: 落ちる
        C('node tools/ledger.mjs', 5000, 0),            // ★対にならず: 落ちる
        C('node tools/mojibake.mjs', 0, 90000),         // ★t0 が無い: 落ちる
        C('git status', 1000, 90000),                   // 対象外
      ];
      const g = tallySolo(cs, [1])[1];
      const row = (id) => g.rows.find((r) => r.id === id);
      t('solo: 単独だけを数える(件数 3)', g.n === 3);
      t('solo: ★混ざった 99 秒を合計に入れない', Math.abs(g.hours - 13 / 3600) < 1e-9);
      t('solo: ★t1 が無い呼び出しを落とす', row('ledger.mjs').n === 0);
      t('solo: ★t0 が無い呼び出しも落とす', row('mojibake.mjs').n === 0);
      t('solo: 道具ごとに割れる', row('check.mjs').n === 2 && row('graph.mjs').n === 1);
      t('solo: 中央値を出す(7s と 3s → 5.0s)', row('check.mjs').med === 5);
      t('solo: 1 件なら中央値はその値', row('graph.mjs').med === 3);
      t('solo: 空なら中央値は null', row('mojibake.mjs').med === null);
      t('solo: 中央値の偶数個は平均', soloMedian([1000, 2000, 3000, 6000]) === 2.5);
      t('solo: 中央値の空は null', soloMedian([]) === null);
      t('solo: ★L0/L1/L2 が同時に出る', Object.keys(tallySolo(cs)).sort().join(',') === '0,1,2');
    }
    // ★L3(素朴)—— M139 の 2.7 時間がどの定義から出たのかを探すために足した
    {
      t('solo: L3 は繋がっていない 1 本だけ当てる', soloOf('node tools/check.mjs --brief', 3) === 'check.mjs');
      t('solo: ★L3 は git と繋がっていたら落とす(L1 は当てる)',
        soloOf('git status && node tools/check.mjs', 3) === null
        && soloOf('git status && node tools/check.mjs', 1) === 'check.mjs');
      t('solo: ★L3 は パイプでも落とす',
        soloOf('node tools/check.mjs | tail -5', 3) === null);
      t('solo: ★L3 は 改行でも落とす',
        soloOf('cd x' + String.fromCharCode(10) + 'node tools/check.mjs', 3) === null);
      t('solo: soloSegments が本数を数える',
        soloSegments('a && b ; c | d') === 4 && soloSegments('a') === 1 && soloSegments('') === 0);
    }
    // ────────────────────────────────────────────────────────
    // ★★M160 —— 保存期間の見張り。★**刈られた形を先に書いてから**道具を試す。
    // ────────────────────────────────────────────────────────
    {
      const O = (at, afterT, usable, first) => ({ at, nRecs: 100, afterT, usable, first, last: 'z' });
      const D0 = '2026-09-08T00:00:00Z';
      const D1 = '2026-09-09T00:00:00Z';
      t('M160: 初回は判定を出さない', m149WatchVerdict(null, O(D0, 2, 1, 'a')).level === 'first');
      t('M160: ★★使える件数が減ったら alarm',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 9, 3, 'a')).level === 'alarm');
      t('M160: ★★T より後が減ったら alarm',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 7, 5, 'a')).level === 'alarm');
      t('M160: ★いちばん古い記録が進んだら warn',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 10, 6, 'b')).level === 'warn');
      t('M160: 増えていれば ok',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 12, 8, 'a')).level === 'ok');
      t('M160: ★横ばい(0 件/日)は warn',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 9, 5, 'a')).level === 'warn');
      t('M160: 速さを 件/日 で出す',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 12, 8, 'a')).perDay === 3);
      t('M160: ★届く日数を出す(134 まで あと 126 を 3 件/日)',
        Math.abs(m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 12, 8, 'a')).etaDays - 126 / 3) < 1e-9);
      t('M160: ★増えないなら届かない(Infinity)',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 9, 5, 'a')).etaDays === Infinity);
      t('M160: ★alarm は warn に負けない(減りと古い側が同時でも alarm)',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 9, 3, 'b')).level === 'alarm');
      t('M160: 観測を組み立てる(窓の端を拾う)', (() => {
        const o = m149Observation([{ ts: '2026-09-01' }, { ts: '2026-09-09' }, { ts: '2026-08-30' }], 7, D0);
        return o.first === '2026-08-30' && o.last === '2026-09-09' && o.nRecs === 3 && o.usable === 7;
      })());
      t('M160: ★T より後だけを afterT に数える', (() => {
        const o = m149Observation([{ ts: '2026-09-06T00:00:00Z' }, { ts: '2026-09-08T00:00:00Z' }], 0, D0);
        return o.afterT === 1;
      })());
      t('M160: 記録が空でも落ちない', (() => {
        const o = m149Observation([], 0, D0);
        return o.first === null && o.last === null && o.afterT === 0;
      })());
      t('M160: 履歴の置き場所は整列で運ばれる ResearchPaper/',
        M149_WATCH_REL.startsWith('ResearchPaper/'));

      // ★★M171(メタ第 34 回)—— 最小間隔の守り。★M170 が踏んだ「14 分で 101 件/日」を再現して塞ぐ。
      const D14 = '2026-09-08T00:14:12Z';                    // ★M170 の実測(14.2 分)
      const D6h = '2026-09-08T06:00:00Z';                    // ★ちょうど 6 時間
      const D6hm = '2026-09-08T05:59:00Z';                   // ★6 時間の 1 分手前
      t('M171: ★★14 分の間隔では perDay を出さない',
        m149WatchVerdict(O(D0, 9, 1, 'a'), O(D14, 11, 2, 'a')).perDay === null);
      t('M171: ★★14 分の間隔では etaDays を出さない',
        m149WatchVerdict(O(D0, 9, 1, 'a'), O(D14, 11, 2, 'a')).etaDays === null);
      t('M171: ★短すぎると言葉で告げる',
        m149WatchVerdict(O(D0, 9, 1, 'a'), O(D14, 11, 2, 'a')).lines.some((l) => l.includes('間隔が短すぎる')));
      t('M171: ★★間隔が短くても level は変わらない(減りは alarm のまま)',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D14, 9, 3, 'a')).level === 'alarm');
      t('M171: ★間隔が短くても古い側が進めば warn',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D14, 10, 6, 'b')).level === 'warn');
      t('M171: ★★6 時間の 1 分手前は まだ出さない',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D6hm, 12, 8, 'a')).perDay === null);
      t('M171: ★ちょうど 6 時間なら出す',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D6h, 12, 8, 'a')).perDay === 12);
      t('M171: ★1 日あけた既定の道は変わらない(3 件/日)',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D1, 12, 8, 'a')).perDay === 3);
      t('M171: ★★短い間隔では「横ばい ⇒ warn」も出さない(0 件/日 と言えないから)',
        m149WatchVerdict(O(D0, 9, 5, 'a'), O(D14, 9, 5, 'a')).level === 'ok');
      t('M171: 最小間隔は 6 時間', M149_WATCH_MIN_DAYS === 0.25);
      t('M171: ★★短い間隔でも「単調に増えている」を押し出さない(第 34 回が実データで踏んだ)',
        m149WatchVerdict(O(D0, 9, 1, 'a'), O(D14, 11, 2, 'a')).lines.some((l) => l.includes('単調に増えている')));
      t('M171: ★見張りの結論が先、間隔の断りが後ろ', (() => {
        const ls = m149WatchVerdict(O(D0, 9, 1, 'a'), O(D14, 11, 2, 'a')).lines;
        return ls.length === 2 && ls[0].includes('単調に増えている') && ls[1].includes('間隔が短すぎる');
      })());
      t('M171: ★alarm のときも間隔の断りは後ろに付く', (() => {
        const v = m149WatchVerdict(O(D0, 9, 5, 'a'), O(D14, 9, 3, 'a'));
        return v.level === 'alarm' && v.lines[v.lines.length - 1].includes('間隔が短すぎる');
      })());

      // ★★M173(メタ第 34 回)—— 「MCP の名指し」の見張り。★分子は第 33 回の検出器のまま。
      t('M173: ★0 件で 10.58% を棄却するのに要る回数は 27',
        binomNeedZero(11 / 104) === 27);
      t('M173: ★率が高いほど早く言える(50% なら 5 回)', binomNeedZero(0.5) === 5);
      t('M173: ★率 0 なら判定不能(null)', binomNeedZero(0) === null);
      t('M173: ★率 1 以上なら判定不能(null)', binomNeedZero(1) === null);
      t('M173: ★★変更後に 1 件でも出たら「まだ起きている」',
        mcpWatchVerdict({ beforeChecks: 104, beforeMiss: 11, afterChecks: 50, afterMiss: 1 }).level === 'まだ起きている');
      t('M173: ★★1 件出たら「言える」には絶対にならない(回数が足りていても)',
        mcpWatchVerdict({ beforeChecks: 104, beforeMiss: 11, afterChecks: 999, afterMiss: 1 }).level !== '言える');
      t('M173: ★0 件のまま必要回数に届けば「言える」',
        mcpWatchVerdict({ beforeChecks: 104, beforeMiss: 11, afterChecks: 27, afterMiss: 0 }).level === '言える');
      t('M173: ★1 回足りなければ「まだ言えない」', (() => {
        const v = mcpWatchVerdict({ beforeChecks: 104, beforeMiss: 11, afterChecks: 26, afterMiss: 0 });
        return v.level === 'まだ言えない' && v.short === 1;
      })());
      t('M173: ★実データの形(照合 1 / 0 件)は あと 26 回', (() => {
        const v = mcpWatchVerdict({ beforeChecks: 104, beforeMiss: 11, afterChecks: 1, afterMiss: 0 });
        return v.level === 'まだ言えない' && v.short === 26;
      })());
      t('M173: ★変更前が空なら基準率が作れない',
        mcpWatchVerdict({ beforeChecks: 0, beforeMiss: 0, afterChecks: 9, afterMiss: 0 }).p0 === null);
      // ★★M192b(メタ第 37 回・持ち場 4)—— 「このままでは溜まらない」を口が自分で言う
      t('M192b: ★変更後 0 件なら stalled(速さ 0)', (() => {
        const a = mcpAccrual({ cut: '2026-09-07T18:16:13Z', lastTs: '2026-09-08T18:16:13Z', afterChecks: 0, short: 26 });
        return a.level === 'stalled' && a.perDay === 0;
      })());
      // ★★★少ない標本から速さを外挿しない —— これが M171 / M174 の守りの再演
      t('M192b: ★★実データの形(変更後 1 回)は stalled で、★速さを出さない', (() => {
        const a = mcpAccrual({ cut: '2026-09-07T18:16:13Z', lastTs: '2026-09-07T20:25:21Z', lastCheckTs: '2026-09-07T18:26:52Z', afterChecks: 1, short: 26 });
        return a.level === 'stalled' && a.perDay === null && a.needDays === null;
      })());
      t('M192b: ★★その場合でも「最後の照合から何時間」は言う', (() => {
        const a = mcpAccrual({ cut: '2026-09-07T18:16:13Z', lastTs: '2026-09-07T20:25:21Z', lastCheckTs: '2026-09-07T18:26:52Z', afterChecks: 1, short: 26 });
        return a.droughtH > 1.9 && a.droughtH < 2.1;
      })());
      t('M192b: ★★「このままでは溜まらない」を必ず口に出す', (() => {
        const a = mcpAccrual({ cut: '2026-09-07T18:16:13Z', lastTs: '2026-09-07T20:25:21Z', afterChecks: 1, short: 26 });
        return a.lines.some((l) => l.includes('このままでは溜まらない'));
      })());
      t('M192b: ★標本が MIN_N 以上なら速さを出す', (() => {
        const a = mcpAccrual({ cut: '2026-09-01T00:00:00Z', lastTs: '2026-09-08T00:00:00Z', afterChecks: 70, short: 26 });
        return a.level === 'accruing' && a.needDays < 3 && a.perDay > 9;
      })());
      t('M192b: ★遅すぎれば stalled(60 日超)', (() => {
        const a = mcpAccrual({ cut: '2026-09-01T00:00:00Z', lastTs: '2026-09-08T00:00:00Z', afterChecks: 5, short: 260 });
        return a.level === 'stalled' && a.needDays > MCP_ACCRUAL_STALL_DAYS;
      })());
      t('M192b: ★もう届いていれば速さを言わない', (() => {
        const a = mcpAccrual({ cut: '2026-09-01T00:00:00Z', lastTs: '2026-09-08T00:00:00Z', afterChecks: 30, short: 0 });
        return a.level === 'reached' && a.lines.length === 0;
      })());
      t('M192b: ★時刻が壊れていても落ちない(stalled として扱う)', (() => {
        const a = mcpAccrual({ cut: 'x', lastTs: null, afterChecks: 3, short: 5 });
        return a.level === 'stalled' && a.days === 0;
      })());
      t('M173: 規約が変わった時刻は 2026-09-07T18:16:13.331Z',
        MCP_WATCH_CUTOFF === '2026-09-07T18:16:13.331Z');
      t('M173: ★分母は「自分で lean_start を呼んだ agent」だけ', (() => {
        const mk = (arr) => arr.map((o) => JSON.stringify(o)).join('\n');
        // 自分では start していない agent(status だけ)は数えない
        const noStart = mk([{ timestamp: 'T1', message: { content: [
          { type: 'tool_use', id: 'u1', name: 'mcp__lean__lean_status', input: {} }] } },
          { timestamp: 'T1', message: { content: [
            { type: 'tool_result', tool_use_id: 'u1', content: 'imports: A, B' }] } }]);
        const withStart = mk([{ timestamp: 'T2', message: { content: [
          { type: 'tool_use', id: 'u0', name: 'mcp__lean__lean_start', input: { imports: ['A'] } }] } },
          { timestamp: 'T2', message: { content: [
            { type: 'tool_result', tool_use_id: 'u0', content: 'imports: A' }] } }]);
        const rows = [{ toolUseId: 'x', description: 'x', s: 0, e: 1 }, { toolUseId: 'y', description: 'y', s: 0, e: 1 }];
        const texts = new Map([['x', noStart], ['y', withStart]]);
        const ev = envCheckEvents(rows, texts);
        return ev.length === 1 && ev[0].id === 'y';
      })());
      t('M173: ★imports を含まない返答は照合に数えない', (() => {
        const mk = (arr) => arr.map((o) => JSON.stringify(o)).join('\n');
        const text = mk([{ timestamp: 'T', message: { content: [
          { type: 'tool_use', id: 'u0', name: 'mcp__lean__lean_start', input: { imports: ['A'] } }] } },
          { timestamp: 'T', message: { content: [
            { type: 'tool_result', tool_use_id: 'u0', content: 'ok(no imports here)' }] } }]);
        return envCheckEvents([{ toolUseId: 'z', description: 'z', s: 0, e: 1 }], new Map([['z', text]])).length === 0;
      })());
    }

    // ── ★★★M166(同時実行数)—— 露出と転帰の純関数 ──────────────────
    {
      const ov = (iv) => overlapStats(iv);
      t('M166: 1 本きりなら atStart=1 / max=1',
        (() => { const m = ov([{ id: 'a', s: 0, e: 100 }]); return m.get('a').atStart === 1 && m.get('a').max === 1; })());
      t('M166: ★重ならない 2 本は互いに 1',
        (() => { const m = ov([{ id: 'a', s: 0, e: 10 }, { id: 'b', s: 20, e: 30 }]);
                 return m.get('a').atStart === 1 && m.get('b').atStart === 1; })());
      t('M166: ★後から始まった側の atStart が 2(先に始まった側は 1)',
        (() => { const m = ov([{ id: 'a', s: 0, e: 100 }, { id: 'b', s: 50, e: 150 }]);
                 return m.get('a').atStart === 1 && m.get('b').atStart === 2; })());
      t('M166: ★★atStart は自分の duration で動かない(同じ開始・違う長さ)',
        (() => { const m = ov([{ id: 'a', s: 0, e: 1000 }, { id: 'b', s: 0, e: 10 }]);
                 return m.get('a').atStart === m.get('b').atStart; })());
      t('M166: ★max は後から重なった本を拾う(atStart=1 でも max=2)',
        (() => { const m = ov([{ id: 'a', s: 0, e: 100 }, { id: 'b', s: 50, e: 150 }]);
                 return m.get('a').atStart === 1 && m.get('a').max === 2; })());
      t('M166: ★mean は時間平均(半分だけ 2 本 ⇒ 1.5)',
        (() => { const m = ov([{ id: 'a', s: 0, e: 100 }, { id: 'b', s: 50, e: 150 }]);
                 return Math.abs(m.get('a').mean - 1.5) < 1e-9; })());
      t('M166: ★3 本が重なれば 3',
        (() => { const m = ov([{ id: 'a', s: 0, e: 100 }, { id: 'b', s: 1, e: 100 }, { id: 'c', s: 2, e: 100 }]);
                 return m.get('c').atStart === 3 && m.get('a').max === 3; })());
      t('M166: 区画は 4 以上を 1 つにまとめる',
        concBin(1) === 1 && concBin(3) === 3 && concBin(4) === 4 && concBin(9) === 4);
      t('M166: 占有 —— 1 本きりなら k=1 に全部',
        (() => { const o = occupancy([{ s: 0, e: 100 }]); return o.byK.get(1) === 100 && o.total === 100; })());
      t('M166: 占有 —— 半分重なれば k=1 が 100 / k=2 が 50', (() => {
        const o = occupancy([{ s: 0, e: 100 }, { s: 50, e: 150 }]);
        return o.byK.get(1) === 100 && o.byK.get(2) === 50 && o.total === 150;
      })());
      t('M166: 占有 —— ★隙間は数えない(誰も走っていない時間)', (() => {
        const o = occupancy([{ s: 0, e: 10 }, { s: 90, e: 100 }]);
        return o.total === 20 && o.span === 100 && !o.byK.has(0);
      })());
      t('M166: 占有 —— ★端が接するだけなら重ならない', (() => {
        const o = occupancy([{ s: 0, e: 50 }, { s: 50, e: 100 }]);
        return !o.byK.has(2) && o.byK.get(1) === 100;
      })());
      t('M166: 占有 —— 3 本の入れ子', (() => {
        const o = occupancy([{ s: 0, e: 30 }, { s: 10, e: 30 }, { s: 20, e: 30 }]);
        return o.byK.get(1) === 10 && o.byK.get(2) === 10 && o.byK.get(3) === 10;
      })());

      // 本文の組み立て(★実物と同じ形。★試験のために最小限)
      const L = (o) => JSON.stringify(o);
      const use = (id, name, input, ts) => L({ timestamp: ts, message: { role: 'assistant', content: [{ type: 'tool_use', id, name, input }] } });
      const res = (id, body, ts) => L({ timestamp: ts, message: { role: 'user', content: [{ type: 'tool_result', tool_use_id: id, content: body }] } });
      const T0 = '2026-09-08T00:00:00.000Z', T1 = '2026-09-08T00:00:30.000Z';

      t('M166(a): ★配管の失敗を拾う(MCP error)',
        concurrencyOutcomesFromText([
          use('1', 'mcp__abc3-lean__lean_check', { code: 'x' }, T0),
          res('1', 'MCP error -32000: no running lean session', T1),
        ].join('\n')).mcpInfra === 1);
      // ★★★感度の試験。★これが無いと「0 件」が嘘になる(メタ第 33 回が実際に踏んだ)。
      for (const f of MCP_INFRA_FIXTURES) {
        t(`M166(a) 感度: ${f.want ? '鳴る' : '黙る'} 「${f.s.slice(0, 26).replace(/\s+/g, ' ')}…」`,
          MCP_INFRA_RE.test(f.s) === f.want);
      }
      t('M166(a): ★★事前登録した v1 は「REPL は処理中」に鳴らない(★この失敗を記録として残す)',
        MCP_INFRA_RE_V1.test('エラー: REPL は処理中(直列にしか使えない)') === false);
      t('M166(a): ★取り合いだけを別に数える(kind.busy)',
        concurrencyOutcomesFromText([use('1', 'mcp__abc3-lean__lean_check', {}, T0),
          res('1', 'エラー: REPL は処理中(直列にしか使えない)', T1)].join('\n')).kind.busy === 1);
      t('M166(a): ★時間切れは取り合いに数えない(kind.timeout)', (() => {
        const o = concurrencyOutcomesFromText([use('1', 'mcp__abc3-lean__lean_check', {}, T0),
          res('1', 'エラー: 600 秒で応答が無いので REPL を落とした。', T1)].join('\n'));
        return o.kind.busy === 0 && o.kind.timeout === 1 && o.mcpInfra === 1;
      })());
      t('M166(a): ★lean_reset の成功は配管の失敗ではない',
        concurrencyOutcomesFromText([use('1', 'mcp__abc3-lean__lean_reset', {}, T0),
          res('1', '[{"type":"text","text":"REPL を落とした(再生用の控えも捨てた)。"}]', T1)].join('\n')).mcpInfra === 0);
      t('M166(a): ★★Lean の型エラーは配管に数えない(ここを混ぜると全部 1 になる)',
        concurrencyOutcomesFromText([
          use('1', 'mcp__abc3-lean__lean_check', { code: 'x' }, T0),
          res('1', 'error: unknown identifier "foo"\nエラー 1 件', T1),
        ].join('\n')).mcpInfra === 0);
      t('M166(a): ★MCP でない道具の失敗も数えない',
        concurrencyOutcomesFromText([
          use('1', 'Bash', { command: 'ls' }, T0),
          res('1', 'MCP error: connection closed', T1),
        ].join('\n')).mcpInfra === 0);
      t('M166(a): lean_start を数える',
        concurrencyOutcomesFromText(use('1', 'mcp__abc3-lean__lean_start', {}, T0)).leanStarts === 1);
      t('M166(c): ★lake build の実時間を tool_use → tool_result の時刻差で測る',
        (() => { const o = concurrencyOutcomesFromText([
          use('1', 'Bash', { command: 'lake build ABC3' }, T0), res('1', 'ok', T1)].join('\n'));
          return o.lakeCalls === 1 && o.lakeWaitMs === 30000; })());
      t('M166(c): ★build.mjs も lake の待ちに数える',
        concurrencyOutcomesFromText(use('1', 'Bash', { command: 'node tools/build.mjs ABC3' }, T0)).lakeCalls === 1);
      t('M166(c): ★lake でない Bash は数えない',
        concurrencyOutcomesFromText(use('1', 'Bash', { command: 'node tools/check.mjs --brief' }, T0)).lakeCalls === 0);
      t('M166(b): ★Write した .lean を拾う(Read は拾わない)',
        (() => { const o = concurrencyOutcomesFromText([
          use('1', 'Write', { file_path: 'D:/x/A.lean', content: '' }, T0),
          use('2', 'Read', { file_path: 'D:/x/B.lean' }, T0)].join('\n'));
          return o.wrote.size === 1 && o.wrote.has('d:/x/a.lean'); })());
      t('M166(b): ★lean-idioms.md も衝突の対象に入れる',
        concurrencyOutcomesFromText(use('1', 'Edit', { file_path: 'D:/x/tools/lean-idioms.md' }, T0)).wrote.size === 1);
      t('M166(b): ★★重なっていなければ同じ本を書いても数えない', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:01:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' },
                   { toolUseId: 'b', ts: '2026-09-08T00:10:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' }];
        const tx = new Map([['a', use('1', 'Write', { file_path: 'X.lean' }, T0)],
                            ['b', use('2', 'Write', { file_path: 'X.lean' }, T0)]]);
        return concurrencyRows(R, tx).every((r) => r.sharedFile === 0);
      })());
      t('M166(b): ★★重なって同じ本を書いたら両方 1', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:01:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' },
                   { toolUseId: 'b', ts: '2026-09-08T00:01:30Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' }];
        const tx = new Map([['a', use('1', 'Write', { file_path: 'X.lean' }, T0)],
                            ['b', use('2', 'Write', { file_path: 'X.lean' }, T0)]]);
        return concurrencyRows(R, tx).every((r) => r.sharedFile === 1);
      })());
      t('M166(b): ★重なっても違う本なら 0', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:01:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' },
                   { toolUseId: 'b', ts: '2026-09-08T00:01:30Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' }];
        const tx = new Map([['a', use('1', 'Write', { file_path: 'X.lean' }, T0)],
                            ['b', use('2', 'Write', { file_path: 'Y.lean' }, T0)]]);
        return concurrencyRows(R, tx).every((r) => r.sharedFile === 0);
      })());
      t('M166: ★lean-prover でない agent は impl の水準に入らない(NaN)', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:01:00Z', durationMs: 60000, agentType: 'general', toolUses: 10, day: 'd' }];
        return Number.isNaN(concurrencyRows(R, new Map())[0].implAtStart);
      })());
      t('M166(d): msPerTool = duration / tool_uses', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:01:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 10, day: 'd' }];
        return concurrencyRows(R, new Map())[0].msPerTool === 6000;
      })());
      t('M166: ★族は 4 本(後から欄を足していない)', FAMILY_C.length === 4);
      t('M166: ★族の欄はこの 4 つで固定',
        FAMILY_C.map((f) => f.key).join(',') === 'mcpInfra,sharedFile,lakeWaitMs,msPerTool');
      t('M166: ★件数が MIN_N 未満なら検定しない', (() => {
        const rs = Array.from({ length: 5 }, (_, i) => ({ implAtStart: i % 3 + 1, mcpInfra: i, sharedFile: i, lakeWaitMs: i, msPerTool: i }));
        return concurrencyLadder(rs, 'implAtStart').rows.length === 0;
      })());
      t('M166: ★分散ゼロの欄は検定しない(全部 0 の欄)', (() => {
        const rs = Array.from({ length: 20 }, (_, i) => ({ implAtStart: i % 4 + 1, mcpInfra: 0, sharedFile: i, lakeWaitMs: i, msPerTool: i }));
        const l = concurrencyLadder(rs, 'implAtStart');
        return l.skipped.some((s) => s.key === 'mcpInfra') && l.rows.every((r) => r.key !== 'mcpInfra');
      })());
      t('M166: 事象 —— ★取り合いの瞬間の同時本数を数える', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:02:00Z', durationMs: 120000, agentType: 'lean-prover', toolUses: 5, day: 'd' },
                   { toolUseId: 'b', ts: '2026-09-08T00:02:30Z', durationMs: 120000, agentType: 'lean-prover', toolUses: 5, day: 'd' }];
        const tx = new Map([
          ['a', [use('1', 'mcp__abc3-lean__lean_check', {}, '2026-09-08T00:01:00Z'),
                 res('1', 'エラー: REPL は処理中(直列にしか使えない)', '2026-09-08T00:01:10Z')].join('\n')],
          ['b', use('2', 'mcp__abc3-lean__lean_check', {}, '2026-09-08T00:01:00Z')]]);
        const ev = busyEvents(concurrencyRows(R, tx), tx);
        return ev.length === 1 && ev[0].agents === 2 && ev[0].mcpAgents === 2;
      })());
      t('M166: 事象 —— ★MCP を使わない agent は mcpAgents に数えない', (() => {
        const R = [{ toolUseId: 'a', ts: '2026-09-08T00:02:00Z', durationMs: 120000, agentType: 'lean-prover', toolUses: 5, day: 'd' },
                   { toolUseId: 'b', ts: '2026-09-08T00:02:30Z', durationMs: 120000, agentType: 'lean-prover', toolUses: 5, day: 'd' }];
        const tx = new Map([
          ['a', [use('1', 'mcp__abc3-lean__lean_check', {}, '2026-09-08T00:01:00Z'),
                 res('1', 'エラー: REPL は処理中(直列にしか使えない)', '2026-09-08T00:01:10Z')].join('\n')],
          ['b', use('2', 'Bash', { command: 'ls' }, '2026-09-08T00:01:00Z')]]);
        const ev = busyEvents(concurrencyRows(R, tx), tx);
        return ev.length === 1 && ev[0].agents === 2 && ev[0].mcpAgents === 1;
      })());
      // ── ★★★★★無音のすり替わり(D27 訂正 第 1077 の壊れ方)
      {
        const R2 = (a, b) => [
          { toolUseId: 'a', ts: '2026-09-08T00:05:00Z', durationMs: 300000, agentType: 'lean-prover', toolUses: 5, day: 'd' },
          { toolUseId: 'b', ts: '2026-09-08T00:05:00Z', durationMs: 300000, agentType: 'lean-prover', toolUses: 5, day: 'd' },
        ];
        const start = (id, imports, ts) => use(id, 'mcp__abc3-lean__lean_start', { imports }, ts);
        const status = (id, ts) => use(id, 'mcp__abc3-lean__lean_status', {}, ts);
        const okRes = (id, imports, ts) => res(id,
          `[{"type":"text","text":"起動: あり\\nimports: ${imports.join(', ')}\\n基準環境: 5"}]`, ts);
        const T = (m) => `2026-09-08T00:0${m}:00Z`;

        t('M166 無音: ★一致していれば鳴らない', (() => {
          const tx = new Map([['a', [start('1', ['X'], T(1)), okRes('1', ['X'], T(1)),
            status('2', T(2)), okRes('2', ['X'], T(2))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
        t('M166 無音: ★★★別 agent の imports に差し替わったら鳴る', (() => {
          const tx = new Map([
            ['a', [start('1', ['X'], T(1)), okRes('1', ['X'], T(1)),
                   status('2', T(2)), okRes('2', ['Y', 'Z'], T(2))].join('\n')],
            ['b', [start('3', ['Y', 'Z'], T(1)), okRes('3', ['Y', 'Z'], T(1))].join('\n')]]);
          const ev = envMismatchEvents(concurrencyRows(R2(), tx), tx);
          return ev.length === 1 && ev[0].culprit !== null;
        })());
        t('M166 無音: ★順序が違うだけでは鳴らない(集合で見る)', (() => {
          const tx = new Map([['a', [start('1', ['X', 'Y'], T(1)),
            status('2', T(2)), okRes('2', ['Y', 'X'], T(2))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
        t('M166 無音: ★★自分で lean_start を呼んでいなければ数えない(誤報を出さない側)', (() => {
          const tx = new Map([['a', [status('2', T(2)), okRes('2', ['Q'], T(2))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
        t('M166 無音: ★2 回 lean_start したらどちらの形でも鳴らない', (() => {
          const tx = new Map([['a', [start('1', ['X'], T(1)), start('2', ['Y'], T(2)),
            status('3', T(3)), okRes('3', ['X'], T(3)),
            status('4', T(4)), okRes('4', ['Y'], T(4))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
        t('M166 無音: ★★相手が同時に走っていなければ culprit は null', (() => {
          const R3 = [
            { toolUseId: 'a', ts: '2026-09-08T00:05:00Z', durationMs: 120000, agentType: 'lean-prover', toolUses: 5, day: 'd' },
            { toolUseId: 'b', ts: '2026-09-08T09:00:00Z', durationMs: 60000, agentType: 'lean-prover', toolUses: 5, day: 'd' }];
          const tx = new Map([
            ['a', [start('1', ['X'], T(3)), status('2', T(4)), okRes('2', ['Y'], T(4))].join('\n')],
            ['b', start('3', ['Y'], '2026-09-08T08:59:00Z')]]);
          const ev = envMismatchEvents(concurrencyRows(R3, tx), tx);
          return ev.length === 1 && ev[0].culprit === null && ev[0].anyOwner === true;
        })());
        t('M166 無音: ★★imports が空(区切りだけ)でも鳴らない', (() => {
          const tx = new Map([['a', [start('1', ['X'], T(1)),
            status('2', T(2)),
            res('2', '[{"type":"text","text":"起動: あり\\nimports:  ,  ,\\n基準環境: 5"}]', T(2))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
        t('M166 無音: ★imports 行が無い出力では鳴らない', (() => {
          const tx = new Map([['a', [start('1', ['X'], T(1)),
            status('2', T(2)), res('2', '[{"type":"text","text":"起動: なし"}]', T(2))].join('\n')], ['b', '']]);
          return envMismatchEvents(concurrencyRows(R2(), tx), tx).length === 0;
        })());
      }
      t('M166: ★Holm は m = 使えた欄の数で掛かる', (() => {
        const rs = Array.from({ length: 30 }, (_, i) => ({ implAtStart: i % 4 + 1, mcpInfra: i % 7, sharedFile: i % 5, lakeWaitMs: i % 3, msPerTool: i % 11 }));
        const l = concurrencyLadder(rs, 'implAtStart');
        return l.rows.length === 4 && l.rows.every((r) => r.holm >= r.p - 1e-12);
      })());
    }
  }

  // --- ★M179 supply(メタ第 35 回)
  {
    t('M179: toRel は Windows の絶対パスを rel にする',
      toRel('D:\\Math_ABC3\\lean\\ABC3\\Found\\PGC\\X.lean') === 'Found/PGC/X.lean');
    t('M179: toRel は / のパスも同じ', toRel('lean/ABC3/Skeleton/PGC/Section1.lean') === 'Skeleton/PGC/Section1.lean');
    t('M179: toRel は worktree の中でも rel を取る',
      toRel('/d/Math_ABC3/.claude/worktrees/w1/lean/ABC3/Found/A.lean') === 'Found/A.lean');
    t('M179: toRel は .lean でなければ null', toRel('lean/ABC3/Found/A.md') === null);
    t('M179: toRel は lean/ABC3 の外なら null', toRel('tools/graph.mjs') === null);
    t('M179: toRel は文字列でなければ null', toRel(null) === null);
    const jl = [
      JSON.stringify({ message: { content: [
        { type: 'tool_use', name: 'Write', input: { file_path: 'lean/ABC3/Found/A.lean' } },
        { type: 'tool_use', name: 'Read',  input: { file_path: 'lean/ABC3/Found/B.lean' } },
      ] } }),
      JSON.stringify({ message: { content: [{ type: 'tool_use', name: 'Edit', input: { file_path: 'lean/ABC3/Skeleton/C.lean' } }] } }),
      JSON.stringify({ message: { content: [{ type: 'tool_use', name: 'Bash', input: { command: 'lean/ABC3/Found/D.lean' } }] } }),
      'これは JSON ではない',
    ].join('\n');
    t('M179: Write / Edit した .lean だけを取る',
      JSON.stringify(writtenLeanRels(jl)) === JSON.stringify(['Found/A.lean', 'Skeleton/C.lean']));
    t('M179: Read は数えない', !writtenLeanRels(jl).includes('Found/B.lean'));
    t('M179: Bash の中の path は数えない', !writtenLeanRels(jl).includes('Found/D.lean'));
    t('M179: 壊れた行で落ちない', writtenLeanRels('{{{').length === 0);
    const FR = new Map([['S/A.lean', { startable: true }], ['S/B.lean', { startable: false }]]);
    const TR = new Set(['S/A.lean', 'S/B.lean', 'F/C.lean']);
    t('M179: 着手可能', supplyClass('S/A.lean', FR, TR) === '着手可能');
    t('M179: 着手不可', supplyClass('S/B.lean', FR, TR) === '着手不可');
    t('M179: 木にあるが sorry 無し', supplyClass('F/C.lean', FR, TR) === 'sorry なし');
    t('M179: 木に無い', supplyClass('F/NEW.lean', FR, TR) === '木に無い');
    t('M179: 書いていない', supplyClass(null, FR, TR) === '書いていない');
    const NODES = [
      { mod: 'M.F', rel: 'F/C.lean', imports: ['M.S1'] },
      { mod: 'M.S1', rel: 'S/A.lean', imports: ['M.S2'] },
      { mod: 'M.S2', rel: 'S/B.lean', imports: [] },
    ];
    t('M179: 上流の着手可能を辿る',
      JSON.stringify(startableUpstream('F/C.lean', NODES, new Set(['S/A.lean']))) === JSON.stringify(['S/A.lean']));
    t('M179: 自分自身は返さない',
      JSON.stringify(startableUpstream('S/A.lean', NODES, new Set(['S/A.lean']))) === JSON.stringify([]));
    t('M179: 2 段上でも辿る',
      JSON.stringify(startableUpstream('F/C.lean', NODES, new Set(['S/B.lean']))) === JSON.stringify(['S/B.lean']));
    t('M179: 木に無い rel は空', startableUpstream('X.lean', NODES, new Set(['S/A.lean'])).length === 0);
    // ★★2026-09-08: 「自分自身を返さない」の見張りは**環が無いと効かない**(突然変異 S7 が素通りした)。
    t('M179: ★環の中でも自分自身は返さない', startableUpstream('S/A.lean',
      [{ mod: 'M.F', rel: 'F/C.lean', imports: ['M.S1'] }, { mod: 'M.S1', rel: 'S/A.lean', imports: ['M.F'] }],
      new Set(['S/A.lean'])).length === 0);
    t('M179: 環があっても止まる', startableUpstream('F/C.lean',
      [{ mod: 'M.F', rel: 'F/C.lean', imports: ['M.S1'] }, { mod: 'M.S1', rel: 'S/A.lean', imports: ['M.F'] }],
      new Set(['S/A.lean'])).length === 1);
    const TL = supplyTally([
      { klass: '着手可能', upstream: ['S/Z.lean'] },   // ★★一覧に載っていた本の上流は数えない(S9)
      { klass: '木に無い', upstream: ['S/A.lean'] },
      { klass: '木に無い', upstream: ['S/A.lean'] },
      { klass: 'sorry なし', upstream: ['S/A.lean', 'S/B.lean'] },
    ]);
    t('M179: 配った本数', TL.n === 4);
    t('M179: 一覧に載っていた本数', TL.listed === 1);
    t('M179: 割合', close(TL.listedPct, 25));
    t('M179: 上流の集計は一覧の外だけを数える', TL.byAncestor[0][0] === 'S/A.lean' && TL.byAncestor[0][1] === 3);
    t('M179: ★一覧に載っていた本の上流は 1 件も混ぜない', !TL.byAncestor.some(([k]) => k === 'S/Z.lean'));
    t('M179: 分類は 5 つ', SUPPLY_CLASSES.length === 5);
    t('M179: import を読む',
      JSON.stringify(parseImports('import ABC3.A\nimport ABC3.B\n\ntheorem x := 1')) === JSON.stringify(['ABC3.A', 'ABC3.B']));
    t('M179: import の後ろに語があれば取らない', parseImports('import ABC3.A -- なにか').length === 0);
    t('M179: import で始まらない行は取らない', parseImports('  -- import ABC3.A').length === 0);
    t('M179: brief は最初の user の文字列', briefTextOf([
      JSON.stringify({ message: { role: 'assistant', content: 'x' } }),
      JSON.stringify({ message: { role: 'user', content: 'これが brief。Skeleton/PGC/Section1.lean' } }),
      JSON.stringify({ message: { role: 'user', content: '2 つ目' } }),
    ].join('\n')).startsWith('これが brief'));
    t('M179: user が無ければ空', briefTextOf('{}') === '');
    t('M179: 字面から Skeleton を拾う',
      JSON.stringify(mentionedSkeletons('… Skeleton/PGC/Section1.lean と Skeleton\\PGC\\Section2.lean …'))
      === JSON.stringify(['Skeleton/PGC/Section1.lean', 'Skeleton/PGC/Section2.lean']));
    t('M179: 重複は 1 つに畳む', mentionedSkeletons('Skeleton/A/B.lean Skeleton/A/B.lean').length === 1);
    t('M179: Found は拾わない', mentionedSkeletons('Found/PGC/X.lean').length === 0);
    const TL2 = supplyTally([
      { klass: '木に無い', upstream: [], briefNodes: ['S/A.lean'], bodyNodes: ['S/A.lean', 'S/B.lean'] },
      { klass: '木に無い', upstream: [], briefNodes: [], bodyNodes: ['S/A.lean'] },
    ]);
    t('M179: brief の水路', TL2.briefHit === 1 && TL2.byBrief[0][1] === 1);
    t('M179: 本文の水路', TL2.bodyHit === 2 && TL2.byBody[0][0] === 'S/A.lean' && TL2.byBody[0][1] === 2);
    t('M179: 合成節点は rel を持ち、import を辿れる',
      JSON.stringify(startableUpstream('F/NEW.lean',
        [synthNode('F/NEW.lean', 'import M.S1'), { mod: 'M.S1', rel: 'S/A.lean', imports: [] }],
        new Set(['S/A.lean']))) === JSON.stringify(['S/A.lean']));
  }

  console.log(`\nselftest: ${ok}/${ok + ng}`);
  return ng === 0;
}

// ════════════════════════════════════════════════════════════════════
// 5b. 診断の族 × 費用(`--diag`)—— ★★事前登録(メタ第 27 回)
// ════════════════════════════════════════════════════════════════════
//
// ★何を測るか / ★★何は測れないか(先に書く)
// ------------------------------------------------
// ★測れない: **「族 → 時間」の因果**。難しいノードは長くもあり診断も多い(交絡)。
//   ★測れるのはせいぜい「族の混合比が違う agent の間で、**1 往復あたりの時間**が違うか」まで。
//   ⇒ 出す判定は「言える / 言えない」だけ。★向き(どちらが高い)は断定しない。
// ★★さらに強い制約(メタ第 27 回の実測): `tree` の診断 4,379 件のうち **4,176 件(95%)は
//   本体セッションが出したもの**で、★**本体には duration_ms が無い**(agent ではないため)。
//   ⇒ 費用と突き合わせられるのは残り 203 件だけ。★これは道具の不備ではなく**観測点の構造**である。
//
// ★事前登録(★データを見る前にここへ焼く。後から欄を足さない)
// ------------------------------------------------
//   母集団 : `agentType === 'lean-prover'`(= 実装者)かつ `toolUses >= 1`。
//            ★役割で切る理由: math-planner / Explore は Lean を叩かないので
//              「診断ゼロ かつ 速い」が構造的に決まり、族の効果と役割が完全に交絡する。
//   結果 Y : `duration_ms / tool_uses`(1 往復あたりのミリ秒)。★中央値で 2 値化(high = Y > 中央値)。
//   曝露 X : その agent が族 f の診断を **1 件以上**出したか。
//   族     : DIAG_FAMILIES(5 つ。★M115 の tree 上位 5 族をそのまま固定)。
//   スコープ: A = `tree` のみ(主) / B = `tree` + `scratch`(副)。★`doc` と `meta` は常に外す。
//   検定   : Fisher 正確検定(両側)。★多重比較は **m = 5 族 × 2 スコープ = 10** の Bonferroni。
//   最小件数: 曝露あり群・曝露なし群ともに MIN_N(8)件以上。★未満なら**検定しない**(M112 と同じ作法)。
export const DIAG_FAMILIES = [
  { id: 'D1', key: 'failed to synthesize' },
  { id: 'D2', key: 'unsolved goals' },
  { id: 'D3', key: 'Did not find an occurrence of the pattern' },
  { id: 'D4', key: 'Application type mismatch' },
  { id: 'D5', key: 'Type mismatch' },
];
export const DIAG_SCOPES = [
  { id: 'A', label: 'tree のみ(木の本を検査して出た失敗)', srcs: ['tree'] },
  { id: 'B', label: 'tree + scratch(REPL の断片も入れる)', srcs: ['tree', 'scratch'] },
];
export const DIAG_TYPE = 'lean-prover';
export const DIAG_ALPHA = 0.05;
export const DIAG_M = DIAG_FAMILIES.length * DIAG_SCOPES.length;   // = 10

/** ★対数ガンマ(Lanczos)。Fisher を n が大きくても引けるようにするため。 */
const LGC = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
function lgamma(z) {
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  z -= 1; let x = 0.99999999999980993;
  for (let i = 0; i < 8; i++) x += LGC[i] / (z + i + 1);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
const lfact = (n) => lgamma(n + 1);
/**
 * Fisher の正確検定(両側・確率法)。★`unverified.mjs` の `fisherLog` と同じ式。
 * ★あちらは厳密版(BigInt)と食い違わないことを selftest で確かめてある。ここでは
 *   ★既知の値(お茶の実験 [[3,1],[1,3]] = 0.4857 など)で較正する。
 */
export function fisher(a, b, c, d) {
  const r1 = a + b, r2 = c + d, k = a + c, n = a + b + c + d;
  if (!r1 || !r2 || !k || k === n) return 1;
  const base = lfact(r1) + lfact(r2) + lfact(k) + lfact(n - k) - lfact(n);
  const lp = (x) => base - lfact(x) - lfact(r1 - x) - lfact(k - x) - lfact(r2 - k + x);
  const obs = lp(a); let sum = 0;
  for (let x = Math.max(0, k - r2); x <= Math.min(k, r1); x++) {
    const v = lp(x);
    if (v <= obs + 1e-9) sum += Math.exp(v);
  }
  return Math.min(1, sum);
}

/**
 * ★「あと何件あれば言えるか」—— 観測した**割合を固定**したまま、均衡した群サイズを増やして
 * Fisher p が α' を割る最小の 1 群あたり件数を返す。★多重比較込みの α' を渡すこと。
 * ★これは検出力計算ではない(割合が本当にそうなら、という条件付きの最小 n)。★見つからなければ null。
 */
export function needBalancedN(a, b, c, d, alpha, cap = 4000) {
  const n1 = a + b, n2 = c + d;
  if (!n1 || !n2) return null;
  const p1 = a / n1, p2 = c / n2;
  if (Math.abs(p1 - p2) < 1e-12) return null;      // 割合が同じなら何件増やしても言えない
  for (let n = MIN_N; n <= cap; n++) {
    const A = Math.round(p1 * n), C = Math.round(p2 * n);
    if (fisher(A, n - A, C, n - C) <= alpha) return n;
  }
  return null;
}

/** digest(idiom-recur v3)を読み、agent(toolUseId)ごとの族の件数にする。★純関数。 */
export function tallyDiagnostics(errors, families = DIAG_FAMILIES, scopes = DIAG_SCOPES) {
  const byAgent = new Map();
  for (const e of errors || []) {
    const ag = e.ag || '';
    if (!ag) continue;                              // ★本体セッション(agent ではない)は突き合わせられない
    const src = e.src || 'tree';
    if (src === 'doc' || src === 'meta') continue;  // ★Lean のエラーではない
    const msg = String(e.msg || '');
    let rec = byAgent.get(ag);
    if (!rec) { rec = {}; for (const s of scopes) rec[s.id] = { total: 0, fam: families.map(() => 0) }; byAgent.set(ag, rec); }
    for (const s of scopes) {
      if (!s.srcs.includes(src)) continue;
      rec[s.id].total++;
      families.forEach((f, i) => { if (msg.includes(f.key)) rec[s.id].fam[i]++; });
    }
  }
  return byAgent;
}

/** 中央値で 2 値化する。★high = 値 > 中央値(同値は low 側)。★事前に決めた向き。 */
export function medianSplit(vals) {
  const med = describe(vals).med;
  return { med, isHigh: (v) => v > med };
}

function cmdDiag(recs, digest) {
  const errs = digest?.errors || [];
  const byAgent = tallyDiagnostics(errs);
  const pop = recs.filter(r => r.agentType === DIAG_TYPE && r.toolUses >= 1 && r.durationMs > 0 && r.toolUseId);
  const rows = pop.map(r => ({
    id: r.toolUseId, day: r.day, ms: r.durationMs, tu: r.toolUses,
    y: r.durationMs / r.toolUses,
    diag: byAgent.get(r.toolUseId) || null,
  }));
  const { med, isHigh } = medianSplit(rows.map(r => r.y));
  const dY = rows.length ? describe(rows.map(r => r.y)) : null;

  console.log('== --diag: 診断の族 × 1 往復あたりの費用(★事前登録。向きは断定しない) ==\n');
  console.log(`  母集団      : agentType=${DIAG_TYPE} かつ tool_uses>=1 …… ${rows.length} 件`);
  console.log(`  結果 Y      : duration_ms / tool_uses。中央値 ${(med / 1000).toFixed(1)} 秒/往復(high = これより上)`);
  if (dY) console.log(`  ★Y の散らばり: 四分位 ${(dY.q1 / 1000).toFixed(1)} — ${(dY.q3 / 1000).toFixed(1)} 秒/往復`
    + `(最小 ${(dY.min / 1000).toFixed(1)} / 最大 ${(dY.max / 1000).toFixed(1)})。★下の族ごとの中央値はこの幅と比べて読むこと。`);
  console.log(`  多重比較    : Bonferroni m = ${DIAG_M}(族 ${DIAG_FAMILIES.length} × スコープ ${DIAG_SCOPES.length})  α' = ${(DIAG_ALPHA / DIAG_M).toFixed(4)}`);
  const attributable = errs.filter(e => e.ag && e.src === 'tree').length;
  const treeAll = errs.filter(e => e.src === 'tree').length;
  console.log(`  ★突き合わせられる tree 診断: ${attributable} / ${treeAll} 件`
    + `(残り ${treeAll - attributable} 件は**本体セッション**が出したもので duration_ms が無い)`);

  const P = (n, w) => String(n).padStart(w);
  for (const s of DIAG_SCOPES) {
    console.log(`\n  -- スコープ ${s.id}: ${s.label} --`);
    console.log('   族                              曝露あり  うち高  曝露なし  うち高  Fisher p   補正後   判定');
    for (let i = 0; i < DIAG_FAMILIES.length; i++) {
      const f = DIAG_FAMILIES[i];
      const ex = rows.filter(r => r.diag && r.diag[s.id].fam[i] > 0);
      const un = rows.filter(r => !(r.diag && r.diag[s.id].fam[i] > 0));
      const a = ex.filter(r => isHigh(r.y)).length, b = ex.length - a;
      const c = un.filter(r => isHigh(r.y)).length, d = un.length - c;
      let p = NaN, adj = NaN, verdict;
      if (ex.length < MIN_N || un.length < MIN_N) {
        verdict = `件数不足(要 ${MIN_N})`;
      } else {
        p = fisher(a, b, c, d); adj = Math.min(1, p * DIAG_M);
        verdict = adj <= DIAG_ALPHA ? '★言える' : '言えない';
      }
      console.log(`   ${f.key.slice(0, 30).padEnd(32)}${P(ex.length, 6)}${P(a, 8)}${P(un.length, 9)}${P(c, 8)}`
        + `${isFinite(p) ? P(p.toFixed(4), 10) : P('—', 10)}${isFinite(adj) ? P(adj.toFixed(3), 9) : P('—', 9)}  ${verdict}`);
      // ★記述(★判定ではない。向きの主張でもない)
      const mEx = ex.length ? describe(ex.map(r => r.y)).med / 1000 : NaN;
      const mUn = un.length ? describe(un.map(r => r.y)).med / 1000 : NaN;
      const need = (ex.length >= MIN_N && un.length >= MIN_N) ? needBalancedN(a, b, c, d, DIAG_ALPHA / DIAG_M) : null;
      // ★所要時間と往復数そのものも並べる(★本体の問いは「時間・往復数はどう違うか」だった)。
      //   ★これは記述であって検定ではない。★検定は Y(秒/往復)だけに事前登録してある。
      const mdEx = ex.length ? describe(ex.map(r => r.ms)).med / 60000 : NaN;
      const mdUn = un.length ? describe(un.map(r => r.ms)).med / 60000 : NaN;
      const tuEx = ex.length ? describe(ex.map(r => r.tu)).med : NaN;
      const tuUn = un.length ? describe(un.map(r => r.tu)).med : NaN;
      console.log(`     記述: 秒/往復 中央値 ${isFinite(mEx) ? mEx.toFixed(1) : '—'} 対 ${isFinite(mUn) ? mUn.toFixed(1) : '—'}`
        + `  ・所要 分 ${isFinite(mdEx) ? mdEx.toFixed(1) : '—'} 対 ${isFinite(mdUn) ? mdUn.toFixed(1) : '—'}`
        + `  ・往復 ${isFinite(tuEx) ? tuEx : '—'} 対 ${isFinite(tuUn) ? tuUn : '—'}`
        + `  ・交絡: 診断総数 ${ex.length ? describe(ex.map(r => r.diag[s.id].total)).med : '—'} 対 ${un.length ? describe(un.map(r => (r.diag ? r.diag[s.id].total : 0))).med : '—'}`
        + `  ・★あと何件(1 群): ${need ?? '—'}`);
    }
  }
  // ★標本は独立でない —— 日ごとの束で ICC を出す(M90 / M99 と同じ作法)
  const byDay = new Map();
  for (const r of rows) { const k = r.day || '?'; if (!byDay.has(k)) byDay.set(k, []); byDay.get(k).push(r.y); }
  const ic = icc1([...byDay.values()]);
  console.log(`\n  ★束(日)= ${ic.k} / N = ${ic.N} / ICC ${isFinite(ic.icc) ? ic.icc.toFixed(3) : '—'}`
    + ` / DEFF ${isFinite(ic.deff) ? ic.deff.toFixed(2) : '—'} / 実効 n ${isFinite(ic.nEff) ? ic.nEff.toFixed(1) : '—'}`);
  console.log('  ★★上の「あと何件」は**実効 n ではなく素の n** の話である(束を無視している。実際にはもっと要る)。');
  console.log('  ★★この表から「族 f が時間を食う」とは言えない。★言えるのは「族 f を出した agent と');
  console.log('    出さなかった agent とで、1 往復あたりの時間の高低の分かれ方が違う(違わない)」まで。');
}

// ════════════════════════════════════════════════════════════════════
// 5c. 本体セッションの費用(`--main`)—— ★★事前登録(メタ第 28 回。M125 が残した宿題)
// ════════════════════════════════════════════════════════════════════
//
// ★なぜ要るか(M122 の実測): Lean を叩く仕事の **88%** は本体がやっているのに、
//   本体は agent ではないので `duration_ms` が無い。★だが `tool_use` と `tool_result` は
//   どちらも時刻つきで会話ログに残っており、`tool_use_id` で対にできる(M125 が 36,251 件で確認)。
//
// ★★事前登録(★M125 が「データを見る前にコードへ焼け」と書いたもの。★後から欄を足さない)
// ------------------------------------------------
//   母集団 : 親セッションのログ `~/.claude/projects/<slug>/*.jsonl`(`subagents/` は入れない)。
//            ★`isSidechain: true` の行は外す(子の会話が混ざる場合の保険。M122 は 3 本とも 0 行と実測)。
//   対応   : `tool_use.id` → 後続の `tool_result.tool_use_id`。時刻は**行の** `timestamp`。
//   ★Y3   : 対の経過秒(= 道具そのものの遅さ)。★人の思考・生成の時間は**入らない**。
//   ★Y2   : 塊あたりの呼び出し数(往復数)。
//   束     : ★**ファイル塊** —— 対象が同じ `.lean` である呼び出しを時刻順に並べ、
//            間隔が MAIN_GAP_MIN(30 分)を超えたら別の塊にする。
//            ★「日」でも「診断の前後」でもなく **ファイル** を単位にすると**先に**決めた。理由は
//            本体の問い(「agent に配れたはずの塊はどれか」)の単位が **1 ファイル = 1 持ち場**だから
//            (agent への brief は「このファイルを埋めよ」の形で配られている)。
//   族     : DIAG_FAMILIES(M122 の 5 つ)を**動かさない**(m の意味が濁ると既存の判定が動く)。
//
// ★★★測れないこと(★先に書く。後から言い訳しない)
//   - 「**配れたはず**」は反実仮想であり、ログからは決して測れない。★因果も測れない。
//   - ★代わりに測るのは「**既に配った塊の輪郭に入るか**」だけ。輪郭の閾は実データから取るが、
//     ★規則はここに固定する:
//       (a) 呼び出し数 ≥ 実際に配った `lean-prover` agent の `tool_uses` の第 1 四分位
//       (b) 塊の期間に **他の `.lean` を Write/Edit していない**(= 1 ファイルで閉じている)
//       (c) `lean_check` か `lake build` を 1 回以上含む(= Lean を実際に叩いている)
//     ⇒ ★これは「**似ている**」であって「配れた」ではない。★言えるのはそこまで。
//   - ★Y3 に人の思考時間は入らない。「稼働時間」は隣接呼び出しの間隔を MAIN_ACTIVE_GAP_S(300 秒)で
//     頭打ちにした**代理**であって、実際に机に向かっていた時間ではない。
export const MAIN_GAP_MIN = 30;
export const MAIN_ACTIVE_GAP_S = 300;
export const MAIN_CACHE_V = 2;   // ★v2: heredoc を剥いだ bkS/tgtS を足した(M142)。v1 の cache は捨てられる。

/** `D:\…\lean\ABC3\Found\X.lean` などを `ABC3/Found/X.lean` に正規化する。★純関数。 */
export function normLean(p) {
  const s = String(p || '').replace(/\\/g, '/');
  const i = s.lastIndexOf('ABC3/');
  return i >= 0 ? s.slice(i) : s;
}

const LEAN_IN_TEXT = /[A-Za-z0-9_.\/\\:-]*ABC3[\\/][A-Za-z0-9_.\/\\-]*\.lean/;

/** 呼び出し 1 件の「対象 .lean」。引数に無ければコマンド本文から拾う。無ければ null。 */
export function leanTarget(name, input) {
  const i = input || {};
  for (const k of ['file_path', 'filePath', 'path', 'notebook_path']) {
    const v = i[k];
    if (typeof v === 'string' && /\.lean$/i.test(v)) return normLean(v);
  }
  const txt = String(i.command ?? i.snippet ?? i.code ?? i.pattern ?? '');
  const m = txt.match(LEAN_IN_TEXT);
  return m ? normLean(m[0]) : null;
}

/**
 * ★ヒアドキュメントの本文は「実行した命令」ではない。判定から外す。
 * ★規則は `tools/hooks/bash-guard.mjs` と同じ(あちらは 2026-08-21 の誤爆で入った)。
 * ★★これが要る理由(M136): 本体は `cat > log.md <<'EOF' … lake build ABC3 … EOF` の形で
 *   **記録の本文に自分のコマンドを書き写す**。素の正規表現はそれを呼び出しとして数える。
 * ★★**既定の判定は変えない**(既存の M95/M101/M110/M112/M128 の分母が動くため)。
 *   剥いだほうは `--main` の「剥いだ」欄と `--denominator` にだけ出す。
 */
export function stripHeredocs(src) {
  const lines = String(src ?? '').split('\n');
  const out = [];
  let tag = null;
  for (const line of lines) {
    if (tag !== null) { if (line.trim() === tag) tag = null; continue; }
    out.push(line);
    const m = line.match(/<<-?\s*(?:'([A-Za-z_][A-Za-z0-9_]*)'|"([A-Za-z_][A-Za-z0-9_]*)"|([A-Za-z_][A-Za-z0-9_]*))/);
    if (m) tag = m[1] || m[2] || m[3];
  }
  return out.join('\n');
}

/** Bash の中身の分類(★記述であって検定ではない。★順序が意味を持つ —— 上から先に当てる)。 */
export const BASH_KINDS = [
  { id: 'lake build', re: /\blake\s+build\b/ },
  { id: 'lake env lean', re: /\blake\s+env\s+lean\b/ },
  { id: 'check.mjs', re: /tools[\\/]check\.mjs/ },
  { id: 'graph.mjs', re: /tools[\\/]graph\.mjs/ },
  { id: 'ledger.mjs', re: /tools[\\/]ledger\.mjs/ },
  { id: 'decl-index.mjs', re: /tools[\\/]decl-index\.mjs/ },
  { id: 'leanfile.mjs', re: /tools[\\/]leanfile\.mjs/ },
  { id: 'tools/ その他', re: /tools[\\/][A-Za-z0-9_.-]+\.(mjs|py|js)/ },
  { id: 'git', re: /(^|[\s;|&(])git\s/ },
  { id: 'grep/rg/find', re: /(^|[\s;|&(])(grep|rg|findstr|find)\s/ },
  { id: 'python', re: /(^|[\s;|&(])(python|py311env)/ },
  { id: 'node(その他)', re: /(^|[\s;|&(])node\s/ },
];
export function classifyBash(cmd) {
  const s = String(cmd || '');
  for (const k of BASH_KINDS) if (k.re.test(s)) return k.id;
  return 'その他';
}

// ══════════════════════════════════════════════════════════════════════════
// ★★★M159 —— 「他と混ざらない命令だけ」の**規則を 1 つに決めて焼く**(メタ第 32 回)
// ══════════════════════════════════════════════════════════════════════════
/**
 * ★なぜ要るか: ★**M139 は「他と混ざらない命令だけ」で 2.7 時間と公表したが、
 *   その規則は台帳にもコードにも無く、再現できない。**★第 31 回は自分の定義で数え直して
 *   4.52 時間を得たが、★「2.7 と 4.52 の差には定義の差も混ざる」と正直に書いて終わっている。
 *   ⇒ ★**ここに規則を焼く。**以後は誰でも同じ数字が出せる(`--solo`)。
 *
 * ★★**正典の規則 SOLO-1**(★これを既定とする。理由は下の「なぜ L1 か」)
 * ------------------------------------------------------------------
 *   Bash 呼び出し 1 件が「道具 T の単独実行」であるとは、
 *     (0) `stripHeredocs` で **heredoc の本文を外した**あと(★M146: 記録の本文に
 *         自分のコマンドを書き写す癖があり、素のままだと幽霊を数える)、
 *     (1) `SOLO_SET` の中で **T だけ**が現れ(他の 5 本は現れない)、
 *     (2) `lake build` も `lake env lean` も現れない
 *   の 3 つがすべて成り立つこと。
 *
 * ★**測る量**: その呼び出しの `tool_result` までの実時間(`t1 - t0`)。
 * ★**母集団**: 親セッション(`isSidechain !== true`)の Bash 呼び出しのみ。
 *
 * ★★**なぜ L1 を正典にするか**: ★M153(第 31 回)が使った定義と**同じ**だから。
 *   ★ここで新しい定義を正典にすると、比べる相手(4.52 時間)が無くなる。
 *   ★より厳しい L2 と、より素朴な L0 も**同時に**数えて並べる ⇒
 *   ★「どちらが正しいか」ではなく「**どの定義でいくつか**」が読めるようにする。
 *
 * ★**この規則が数えないもの(既知の穴。隠さず書く)**:
 *   - `&&` で `git` や `grep` が繋がった呼び出しは L1 では**単独**に数える
 *     (それらは `SOLO_SET` に居ないため)。⇒ L2 でも `tools/` 以外は落とさない。
 *   - 同じ道具を 1 命令の中で 2 回叩いた場合も 1 件と数える(★M130 の「2 度建て」と同じ形)。
 *   - `tool_result` が対になっていない呼び出しは母集団から落ちる(時間が測れないため)。
 */
export const SOLO_SET = [
  { id: 'check.mjs', re: /tools[\\/]check\.mjs/ },
  { id: 'mojibake.mjs', re: /tools[\\/]mojibake\.mjs/ },
  { id: 'decl-index.mjs', re: /tools[\\/]decl-index\.mjs/ },
  { id: 'graph.mjs', re: /tools[\\/]graph\.mjs/ },
  { id: 'unverified.mjs', re: /tools[\\/]unverified\.mjs/ },
  { id: 'ledger.mjs', re: /tools[\\/]ledger\.mjs/ },
];
export const SOLO_LAKE = /\blake\s+(?:build|env\s+lean)\b/;
/** ★L2 でだけ使う: 命令文に出てくる `tools/xxx.(mjs|py|js)` の**基底名**を全部拾う。 */
export function soloToolsMentioned(s) {
  const out = new Set();
  const re = /tools[\\/]([A-Za-z0-9_.-]+\.(?:mjs|py|js))/g;
  let m;
  while ((m = re.exec(String(s || '')))) out.add(m[1]);
  return out;
}

/**
 * ★純関数。命令文 1 本が「どの道具の単独実行か」を返す(該当しなければ `null`)。
 * @param {string} cmd   Bash の命令文(素のまま渡してよい。L>=1 は中で剥ぐ)
 * @param {number} level 0=素(剥がない) / 1=★正典 SOLO-1 / 2=厳格(他の tools/ も許さない)
 *                       / 3=素朴(★命令が 1 本だけ。`&&` `;` `|` で繋がっていない)
 */
export function soloOf(cmd, level = 1) {
  const raw = String(cmd ?? '');
  const s = level >= 1 ? stripHeredocs(raw) : raw;
  if (SOLO_LAKE.test(s)) return null;                 // (2) lake は混ぜない
  const hit = SOLO_SET.filter((k) => k.re.test(s));
  if (hit.length !== 1) return null;                  // (1) ちょうど 1 本
  if (level >= 2) {
    const mentioned = soloToolsMentioned(s);
    if (mentioned.size !== 1) return null;            // ★他の tools/ が居たら落とす
    if (/(^|[\s;|&(])(python|py311env)/.test(s)) return null;
  }
  if (level >= 3 && soloSegments(s) !== 1) return null; // ★命令が 1 本だけ
  return hit[0].id;
}

/**
 * ★L3 でだけ使う: 命令文がいくつの「命令」から成るか。
 * ★`&&` `||` `;` `|` と改行で割り、空でない断片を数える。
 * ★これは**近似**である(引用符の中の `;` も割ってしまう)。★近似だと分かるように名前を分けてある。
 */
export function soloSegments(s) {
  return String(s ?? '')
    .split(/&&|\|\||[;|\n]/)
    .map((x) => x.trim())
    .filter(Boolean).length;
}

/** ★中央値(ミリ秒の配列 → 秒)。★空なら null。 */
export function soloMedian(ms) {
  if (!ms.length) return null;
  const a = ms.slice().sort((x, y) => x - y);
  const n = a.length;
  const v = n % 2 ? a[(n - 1) / 2] : (a[n / 2 - 1] + a[n / 2]) / 2;
  return v / 1000;
}

/** ★純関数。呼び出しの配列 → 定義ごと・道具ごとの集計。★I/O をしない ⇒ selftest で較正できる。 */
export function tallySolo(calls, levels = [0, 1, 2]) {
  const out = {};
  for (const lv of levels) {
    const per = {};
    for (const k of SOLO_SET) per[k.id] = [];
    for (const c of calls) {
      if (!c.cmd) continue;
      if (!(c.t0 > 0 && c.t1 > 0 && c.t1 >= c.t0)) continue;
      const id = soloOf(c.cmd, lv);
      if (id) per[id].push(c.t1 - c.t0);
    }
    const rows = SOLO_SET.map((k) => ({
      id: k.id, n: per[k.id].length,
      hours: per[k.id].reduce((a, b) => a + b, 0) / 3600000,
      med: soloMedian(per[k.id]),
    }));
    out[lv] = { rows, n: rows.reduce((a, r) => a + r.n, 0), hours: rows.reduce((a, r) => a + r.hours, 0) };
  }
  return out;
}

/** 親セッションのログを舐めて `tool_use` ↔ `tool_result` の対を作る。★重いので cache を持つ。 */
export async function scanMainCalls(opts = {}) {
  const roots = opts.roots ?? projectRoots();
  const files = [];
  for (const root of roots) {
    if (!fs.existsSync(root)) continue;
    for (const f of fs.readdirSync(root)) if (f.endsWith('.jsonl')) files.push(path.join(root, f));
  }
  files.sort();
  const sig = files.map(f => { const s = fs.statSync(f); return `${f}:${s.size}:${Math.floor(s.mtimeMs)}`; }).join('|');
  const cache = opts.cache ?? null;
  if (cache && fs.existsSync(cache)) {
    try {
      const c = JSON.parse(fs.readFileSync(cache, 'utf8'));
      if (c.v === MAIN_CACHE_V && c.sig === sig) return { ...c, nFiles: files.length, cached: true };
    } catch { /* 作り直す */ }
  }
  const calls = [];
  let uses = 0, sidechain = 0, unpaired = 0;
  for (const f of files) {
    const sess = path.basename(f, '.jsonl').slice(0, 8);
    const pend = new Map();
    const rl = readline.createInterface({ input: fs.createReadStream(f), crlfDelay: Infinity });
    for await (const line of rl) {
      if (!line.includes('"tool_use"') && !line.includes('"tool_result"')) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }
      if (o.isSidechain === true) { sidechain++; continue; }
      const ts = o.timestamp ? Date.parse(o.timestamp) : 0;
      const c = o.message?.content;
      if (!Array.isArray(c)) continue;
      for (const b of c) {
        if (b.type === 'tool_use') {
          uses++;
          const cmd = b.name === 'Bash' ? String(b.input?.command || '') : '';
          // ★`bk`/`tgt` は**素の命令文**(事前登録の分母。既存の判定を動かさない)。
          // ★`bkS`/`tgtS` は**heredoc の本文を外した**もの(M142 の新しい欄。並べて出すだけ)。
          const stripped = cmd ? stripHeredocs(cmd) : '';
          pend.set(b.id, { id: b.id, s: sess, n: String(b.name || '?'), bk: cmd ? classifyBash(cmd) : '',
                           bkS: cmd ? classifyBash(stripped) : '',
                           t0: ts, t1: 0, tgt: leanTarget(b.name, b.input),
                           tgtS: leanTarget(b.name, cmd ? { ...b.input, command: stripped } : b.input) });
        } else if (b.type === 'tool_result') {
          const u = pend.get(b.tool_use_id);
          if (!u) continue;
          pend.delete(b.tool_use_id);
          u.t1 = ts; u.e = b.is_error === true ? 1 : 0;
          calls.push(u);
        }
      }
    }
    unpaired += pend.size;
  }
  calls.sort((a, b) => a.t0 - b.t0);
  const out = { v: MAIN_CACHE_V, sig, uses, sidechain, unpaired, calls };
  if (cache) { fs.mkdirSync(path.dirname(cache), { recursive: true }); fs.writeFileSync(cache, JSON.stringify(out)); }
  return { ...out, nFiles: files.length, cached: false };
}

/**
 * ★M159 —— `--solo` **専用**の走査。★`scanMainCalls` の cache には**手を触れない**。
 *
 * ★★なぜ別の走査にするか(★これは設計上の判断であって手抜きではない):
 *   `scanMainCalls` の記録は `bk`(分類)しか持たず**命令文そのものを捨てている**ので、
 *   規則を後から当て直せない。★命令文を足すには `MAIN_CACHE_V` を上げるしかなく、
 *   ★そうすると **M128 / M139 / M149 が寄りかかっている cache が作り直される**。
 *   台帳は繰り返し「★既定の判定は変えない」と書いている ⇒ ★**別の口・別の cache**にする。
 */
export const SOLO_CACHE_V = 1;
export async function scanSoloCalls(opts = {}) {
  const roots = opts.roots ?? projectRoots();
  const files = [];
  for (const root of roots) {
    if (!fs.existsSync(root)) continue;
    for (const f of fs.readdirSync(root)) if (f.endsWith('.jsonl')) files.push(path.join(root, f));
  }
  files.sort();
  const sig = files.map((f) => { const s = fs.statSync(f); return `${f}:${s.size}:${Math.floor(s.mtimeMs)}`; }).join('|');
  const cache = opts.cache ?? null;
  if (cache && fs.existsSync(cache)) {
    try {
      const c = JSON.parse(fs.readFileSync(cache, 'utf8'));
      if (c.v === SOLO_CACHE_V && c.sig === sig) return { ...c, nFiles: files.length, cached: true };
    } catch { /* 作り直す */ }
  }
  const calls = [];
  let bashUses = 0, sidechain = 0, unpaired = 0;
  for (const f of files) {
    const pend = new Map();
    const rl = readline.createInterface({ input: fs.createReadStream(f), crlfDelay: Infinity });
    for await (const line of rl) {
      if (!line.includes('"tool_use"') && !line.includes('"tool_result"')) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }
      if (o.isSidechain === true) { sidechain++; continue; }
      const ts = o.timestamp ? Date.parse(o.timestamp) : 0;
      const c = o.message?.content;
      if (!Array.isArray(c)) continue;
      for (const b of c) {
        if (b.type === 'tool_use') {
          if (b.name !== 'Bash') continue;            // ★母集団は Bash だけ
          bashUses++;
          // ★規則は走査の時点で当てず、命令文を持って帰る(定義を後から並べ替えられるように)。
          pend.set(b.id, { t0: ts, t1: 0, cmd: String(b.input?.command || '') });
        } else if (b.type === 'tool_result') {
          const u = pend.get(b.tool_use_id);
          if (!u) continue;
          pend.delete(b.tool_use_id);
          u.t1 = ts;
          // ★どの定義にも当たらない命令文は**捨てる**(cache を小さく保つ。L0/L1/L2 のいずれかに当たれば残す)
          if (soloOf(u.cmd, 0) || soloOf(u.cmd, 1) || soloOf(u.cmd, 2)) calls.push(u);
        }
      }
    }
    unpaired += pend.size;
  }
  calls.sort((a, b) => a.t0 - b.t0);
  const out = { v: SOLO_CACHE_V, sig, bashUses, sidechain, unpaired, calls };
  if (cache) { fs.mkdirSync(path.dirname(cache), { recursive: true }); fs.writeFileSync(cache, JSON.stringify(out)); }
  return { ...out, nFiles: files.length, cached: false };
}

/** 呼び出しの種別(★塊の条件 (c) に使う)。 */
export const isLeanCheck = (c) => /lean_check$/.test(c.n);
export const isLakeBuild = (c) => c.bk === 'lake build';
export const isWriteEdit = (c) => c.n === 'Write' || c.n === 'Edit' || c.n === 'MultiEdit' || c.n === 'NotebookEdit';

/** ★ファイル塊(束の単位)。同じ `.lean` の呼び出しを時刻順に並べ、gap 分で切る。★純関数。 */
export function bundleByFile(calls, gapMin = MAIN_GAP_MIN) {
  const byFile = new Map();
  for (const c of calls) {
    if (!c.tgt) continue;
    if (!byFile.has(c.tgt)) byFile.set(c.tgt, []);
    byFile.get(c.tgt).push(c);
  }
  const out = [];
  for (const [file, cs] of byFile) {
    cs.sort((a, b) => a.t0 - b.t0);
    let cur = null;
    for (const c of cs) {
      const end = Math.max(c.t1 || 0, c.t0);
      if (!cur || (c.t0 - cur.end) > gapMin * 60000) { cur = { file, start: c.t0, end, calls: [] }; out.push(cur); }
      cur.calls.push(c);
      cur.end = Math.max(cur.end, end);
    }
  }
  out.sort((a, b) => a.start - b.start);
  return out;
}

/**
 * ★「既に配った塊の輪郭に入るか」。★これは「配れた」の判定ではない(反実仮想は測れない)。
 * ★規則は事前登録(a)(b)(c)。★理由の文字列も返す(落ちた理由が数えられるように)。
 */
export function delegLike(bundle, writeEdits, minCalls) {
  const n = bundle.calls.length;
  if (n < minCalls) return { ok: false, why: `短い(<${minCalls})` };
  if (!bundle.calls.some(c => isLeanCheck(c) || isLakeBuild(c))) return { ok: false, why: 'Lean を叩いていない' };
  let foreign = 0;
  for (const w of writeEdits) {
    if (w.t0 < bundle.start) continue;
    if (w.t0 > bundle.end) break;
    if (w.tgt && w.tgt !== bundle.file) foreign++;
  }
  if (foreign > 0) return { ok: false, why: `他の .lean を ${foreign} 回書いている` };
  return { ok: true, why: '' };
}

/** 稼働時間の代理。★隣接呼び出しの間隔を cap で頭打ちにして足す(★実作業時間ではない)。 */
export function activeMs(calls, capS = MAIN_ACTIVE_GAP_S) {
  const cs = calls.filter(c => c.t0).slice().sort((a, b) => a.t0 - b.t0);
  let sum = 0;
  for (let i = 0; i < cs.length; i++) {
    const end = Math.max(cs[i].t1 || 0, cs[i].t0);
    sum += end - cs[i].t0;
    if (i + 1 < cs.length) { const gap = cs[i + 1].t0 - end; if (gap > 0) sum += Math.min(gap, capS * 1000); }
  }
  return sum;
}

/**
 * digest(idiom-recur v3)の**本体セッション**の診断を、秒で呼び出しに割り当てる。
 * ★同じ秒に複数の呼び出しが終わっているときは対象ファイルで割り、割れなければ **ambiguous** に数える
 *   (★黙って片方に寄せない)。★agent の診断(`ag` あり)はここでは扱わない。
 */
export function attributeDiagnostics(calls, errors, srcs = ['tree'], families = DIAG_FAMILIES) {
  const bySec = new Map();
  for (const c of calls) {
    if (!c.t1) continue;
    const k = Math.floor(c.t1 / 1000);
    if (!bySec.has(k)) bySec.set(k, []);
    bySec.get(k).push(c);
  }
  const byCall = new Map();
  let matched = 0, ambiguous = 0, unmatched = 0;
  for (const e of errors || []) {
    if (e.ag) continue;
    if (!srcs.includes(e.src || 'tree')) continue;
    const cs = bySec.get(e.ts) || [];
    let pick = null;
    if (cs.length === 1) pick = cs[0];
    else if (cs.length > 1) {
      const f = normLean(e.file || '');
      const same = cs.filter(c => c.tgt && f && c.tgt === f);
      if (same.length === 1) pick = same[0];
      else { ambiguous++; continue; }
    }
    if (!pick) { unmatched++; continue; }
    matched++;
    let r = byCall.get(pick.id);
    if (!r) { r = { total: 0, fam: families.map(() => 0) }; byCall.set(pick.id, r); }
    r.total++;
    const msg = String(e.msg || '');
    families.forEach((f2, i) => { if (msg.includes(f2.key)) r.fam[i]++; });
  }
  return { byCall, matched, ambiguous, unmatched };
}

const hr = (ms) => (ms / 3600000).toFixed(1);
const sec1 = (ms) => (ms / 1000).toFixed(2);

async function cmdMain(recs, digest, opts = {}) {
  const t0 = Date.now();
  const sc = await scanMainCalls({ cache: opts.cache });
  const calls = sc.calls.filter(c => c.t1 && c.t1 >= c.t0);
  const dur = calls.map(c => c.t1 - c.t0);
  const D = describe(dur);
  console.log('== --main: 本体セッションの費用(★事前登録。★記述のみ。因果は測れない) ==\n');
  console.log(`  親セッション: ${sc.nFiles} 本 / tool_use ${sc.uses} 件 / 対になった ${calls.length} 件`
    + `(${(100 * calls.length / Math.max(1, sc.uses)).toFixed(1)}%) / 対にならなかった ${sc.unpaired} 件`
    + ` / isSidechain で外した行 ${sc.sidechain}`);
  console.log(`  ★Y3(対の経過秒): 中央値 ${sec1(D.med)} / 四分位 ${sec1(D.q1)} — ${sec1(D.q3)}`
    + ` / 最大 ${sec1(D.max)} / ★合計 ${hr(D.sum)} 時間`);
  const act = activeMs(calls);
  console.log(`  ★稼働の代理(間隔を ${MAIN_ACTIVE_GAP_S} 秒で頭打ち): ${hr(act)} 時間`
    + ` ⇒ ★道具の待ちが占めるのは ${(100 * D.sum / Math.max(1, act)).toFixed(1)}%(★残りは思考・生成・人の待ち。★内訳は測れない)`);
  console.log(`  ★取得 ${((Date.now() - t0) / 1000).toFixed(1)} 秒(${sc.cached ? 'cache' : '走査'})`);

  // -- 1. 道具ごと
  const byTool = new Map();
  for (const c of calls) {
    let r = byTool.get(c.n); if (!r) { r = { n: 0, ms: [], err: 0 }; byTool.set(c.n, r); }
    r.n++; r.ms.push(c.t1 - c.t0); r.err += c.e || 0;
  }
  const rows = [...byTool.entries()].map(([k, r]) => ({ k, n: r.n, d: describe(r.ms), err: r.err }))
    .sort((a, b) => b.d.sum - a.d.sum);
  console.log('\n  -- 1. 道具ごと(★1 回あたりの所要 と 回数。合計の降順) --');
  console.log('   道具                          回数   中央値 秒   四分位 秒        合計 時間   占有%   is_error');
  for (const r of rows) {
    console.log(`   ${padr(r.k.replace(/^mcp__abc3-lean__/, 'mcp:'), 28)}${pad(r.n, 6)}${pad(sec1(r.d.med), 10)}`
      + `${pad(sec1(r.d.q1) + '—' + sec1(r.d.q3), 16)}${pad(hr(r.d.sum), 12)}${pad((100 * r.d.sum / D.sum).toFixed(1), 8)}${pad(r.err, 10)}`);
  }

  // -- 2. Bash の中身
  const byBash = new Map();
  for (const c of calls) {
    if (c.n !== 'Bash') continue;
    let r = byBash.get(c.bk || 'その他'); if (!r) { r = { n: 0, ms: [] }; byBash.set(c.bk || 'その他', r); }
    r.n++; r.ms.push(c.t1 - c.t0);
  }
  const brows = [...byBash.entries()].map(([k, r]) => ({ k, n: r.n, d: describe(r.ms) })).sort((a, b) => b.d.sum - a.d.sum);
  // ★★M142 の新しい欄: heredoc の本文を外したときの回数と合計。★既定の列は動かさない。
  const bashCalls = calls.filter(c => c.n === 'Bash');
  const strN = new Map(), strMs = new Map(), ghostN = new Map(), ghostMs = new Map();
  let gN = 0, gMs = 0;
  for (const c of bashCalls) {
    const k0 = c.bk || 'その他', k1 = c.bkS || 'その他', d = c.t1 - c.t0;
    strN.set(k1, (strN.get(k1) || 0) + 1); strMs.set(k1, (strMs.get(k1) || 0) + d);
    if (k0 !== k1) { ghostN.set(k0, (ghostN.get(k0) || 0) + 1); ghostMs.set(k0, (ghostMs.get(k0) || 0) + d); gN++; gMs += d; }
  }
  console.log('\n  -- 2. Bash の中身(★分類は BASH_KINDS。上から先に当てる) --');
  console.log('   ★「剥いだ」= heredoc の本文を判定から外したとき(M142)。★既定の判定は左の列のまま。');
  console.log('   種類                          回数   中央値 秒   四分位 秒        合計 時間  ★剥いだ回数  ★剥いだ合計h  ★幽霊件数');
  for (const r of brows) {
    console.log(`   ${padr(r.k, 28)}${pad(r.n, 6)}${pad(sec1(r.d.med), 10)}${pad(sec1(r.d.q1) + '—' + sec1(r.d.q3), 16)}${pad(hr(r.d.sum), 12)}`
      + `${pad(strN.get(r.k) ?? 0, 13)}${pad(hr(strMs.get(r.k) ?? 0), 14)}${pad(ghostN.get(r.k) ?? 0, 11)}`);
  }
  const tgtRaw = calls.filter(c => c.tgt).length, tgtStr = calls.filter(c => c.tgtS).length;
  console.log(`   ★幽霊(素で当たり、剥ぐと当たらない)= ${gN} 件 / ${hr(gMs)} 時間`
    + ` —— ★上の族の順序のせいで、下の族(git・check.mjs 等)は素の列で**少なく**出る。`);
  console.log(`   ★対象 .lean が付いた呼び出し: 素 ${tgtRaw} → 剥ぐと ${tgtStr}(差 ${tgtStr - tgtRaw})`
    + ' —— ★下の「ファイル塊」はこの素の側で作ってある。');

  // -- 3. 失敗形ごと(★スコープは M122 が事前登録した A / B をそのまま使う)
  const errs = digest?.errors || [];
  const callById = new Map(calls.map(c => [c.id, c]));
  let at = null;
  console.log('\n  -- 3. 失敗形ごと(★族は M122 の 5 つ。★足していない) --');
  for (const s of DIAG_SCOPES) {
    const a = attributeDiagnostics(calls, errs, s.srcs);
    if (s.id === 'A') at = a;
    const mine = errs.filter(e => !e.ag && s.srcs.includes(e.src || 'tree')).length;
    console.log(`\n   [スコープ ${s.id}] ${s.label}`);
    console.log(`   本体の診断 ${mine} 件 → 呼び出しに割り当てられた ${a.matched}`
      + ` / 同じ秒で割れない ${a.ambiguous} / 対応する呼び出しが無い ${a.unmatched}`);
    console.log('   族                              診断件数  それを出した呼び出し  中央値 秒  合計 時間');
    for (let i = 0; i < DIAG_FAMILIES.length; i++) {
      const f = DIAG_FAMILIES[i];
      let cnt = 0; const ms = [];
      for (const [id, r] of a.byCall) {
        if (!r.fam[i]) continue;
        cnt += r.fam[i];
        const c = callById.get(id); if (c) ms.push(c.t1 - c.t0);
      }
      const d = ms.length ? describe(ms) : null;
      console.log(`   ${padr(f.key.slice(0, 30), 32)}${pad(cnt, 8)}${pad(ms.length, 20)}`
        + `${pad(d ? sec1(d.med) : '—', 11)}${pad(d ? hr(d.sum) : '—', 10)}`);
    }
  }
  console.log('   ★この「合計 時間」は**その診断が出た呼び出しそのものの所要**であって、');
  console.log('     ★その失敗形に費やした時間ではない(直し方を考える時間も、直す往復も入らない)。');

  // -- 4. ファイル塊 と 「配った塊の輪郭」
  const bundles = bundleByFile(calls);
  const withTgt = calls.filter(c => c.tgt).length;
  const pop = recs.filter(r => r.agentType === DIAG_TYPE && r.toolUses >= 1);
  const minCalls = pop.length ? Math.round(quantile(pop.map(r => r.toolUses), 0.25)) : NaN;
  const writeEdits = calls.filter(c => isWriteEdit(c) && c.tgt).sort((a, b) => a.t0 - b.t0);
  const judged = bundles.map(b => ({ b, v: delegLike(b, writeEdits, minCalls) }));
  const like = judged.filter(x => x.v.ok);
  const bcalls = (b) => b.calls.length;
  const bms = (b) => b.calls.reduce((s, c) => s + (c.t1 - c.t0), 0);
  console.log('\n  -- 4. ファイル塊(★束の単位。gap ' + MAIN_GAP_MIN + ' 分) --');
  console.log(`   対象 .lean が取れた呼び出し ${withTgt} / ${calls.length} 件 ⇒ 塊 ${bundles.length} 本`
    + ` / 触った .lean ${new Set(bundles.map(b => b.file)).size} 本`);
  if (bundles.length) {
    const dB = describe(bundles.map(bcalls));
    console.log(`   ★Y2(塊あたりの往復数): 中央値 ${dB.med} / 四分位 ${dB.q1} — ${dB.q3} / 最大 ${dB.max}`);
  }
  console.log(`   ★輪郭の閾(a): 配った lean-prover ${pop.length} 件の tool_uses の第 1 四分位 = ${minCalls} 往復`);
  console.log(`   ★輪郭に入る塊: ${like.length} / ${bundles.length} 本`
    + `(呼び出し ${like.reduce((s, x) => s + bcalls(x.b), 0)} 件 / 対の合計 ${hr(like.reduce((s, x) => s + bms(x.b), 0))} 時間)`);
  const why = new Map();
  for (const x of judged) if (!x.v.ok) { const k = x.v.why.replace(/\d+/g, 'N'); why.set(k, (why.get(k) || 0) + 1); }
  console.log('   ★落ちた理由: ' + [...why.entries()].sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k} ${v}`).join(' / '));
  const lim = opts.limit ?? 15;
  console.log(`\n   -- 輪郭に入った塊(往復の多い順 ${lim} 本) --`);
  console.log('   開始(UTC)          往復  対の分  幅の分  check  build  診断  ファイル');
  for (const x of like.slice().sort((a, b) => bcalls(b.b) - bcalls(a.b)).slice(0, lim)) {
    const b = x.b;
    const nd = b.calls.reduce((s, c) => s + (at.byCall.get(c.id)?.total || 0), 0);
    console.log(`   ${padr(new Date(b.start).toISOString().replace('T', ' ').slice(0, 16), 18)}`
      + `${pad(b.calls.length, 5)}${pad((bms(b) / 60000).toFixed(1), 8)}${pad(((b.end - b.start) / 60000).toFixed(1), 8)}`
      + `${pad(b.calls.filter(isLeanCheck).length, 7)}${pad(b.calls.filter(isLakeBuild).length, 7)}${pad(nd, 6)}  ${b.file}`);
  }
  console.log('   ★★この表は「**配れたはず**」ではない。「配った塊と同じ形をしている」までである。');
  console.log('     ★判断(本当に配れるか)は brief が書けるかどうかで決まり、それはログからは測れない。');

  // -- 5. lake build の反復
  const lb = calls.filter(isLakeBuild);
  if (lb.length) {
    const d = describe(lb.map(c => c.t1 - c.t0));
    let runs = 0, maxRun = 0, cur = 0;
    for (const c of calls) {
      if (isLakeBuild(c)) { cur++; maxRun = Math.max(maxRun, cur); }
      else { if (cur >= 3) runs++; cur = 0; }
    }
    if (cur >= 3) runs++;
    console.log(`\n  -- 5. lake build の反復 --`);
    console.log(`   ${lb.length} 回 / 中央値 ${sec1(d.med)} 秒 / 合計 ${hr(d.sum)} 時間`
      + ` / ★間に何も挟まず 3 回以上続いた塊 ${runs} 本(最長 ${maxRun} 連)`);
  }
}

// ════════════════════════════════════════════════════════════════════
// 6. CLI
// ════════════════════════════════════════════════════════════════════
async function main() {
  const argv = process.argv.slice(2);
  const has = (f) => argv.includes(f);
  const val = (f, d) => { const i = argv.indexOf(f); return i >= 0 && argv[i + 1] ? argv[i + 1] : d; };

  if (has('--selftest')) { process.exit(selftest() ? 0 : 1); }

  const repoRoot = path.resolve(path.join(new URL('.', import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1'), '..'));
  let recs = collect();
  const type = val('--type', null);
  const day = val('--day', null);
  const name = val('--name', null);
  if (type) recs = recs.filter(r => r.agentType === type);
  if (day) recs = recs.filter(r => r.day === day);
  if (name) recs = recs.filter(r => new RegExp(name).test(r.description));
  recs = filterByTime(recs, val('--since', null), val('--until', null));
  const ms = val('--ms', null);           // ★duration_ms のリストで拾う(過去の表を再現するため)
  if (ms) {
    const want = new Set(ms.split(/[,\s]+/).filter(Boolean).map(Number));
    recs = recs.filter(r => want.has(r.durationMs));
  }

  const wantCov = has('--stats') || has('--explain') || has('--estimate') || has('--cov') || has('--cost');
  if (wantCov) {
    for (const r of recs) Object.assign(r, covariates(r, repoRoot));
  }

  if (has('--json')) { console.log(JSON.stringify(recs, null, 1)); return; }

  const filt = `type=${type ?? '全部'} day=${day ?? '全部'}${name ? ' name=/' + name + '/' : ''}`;
  console.log(`== agent-timing —— 通知 ${recs.length} 件(${filt}) ==\n`);

  if (has('--stats')) { cmdStats(recs); console.log(''); cmdIcc(recs); return; }
  if (has('--explain')) { cmdExplain(recs); cmdIcc(recs); return; }
  if (has('--estimate')) { cmdEstimate(recs); return; }
  if (has('--denominator')) { cmdDenominator(recs, repoRoot); return; }
  // ★M160 は `--m149` より先に見る(`--m149-watch` は `--m149` を含む文字列ではないが、
  //   将来 `has()` が前方一致になっても取り違えないように順序で守る)。
  if (has('--supply')) { cmdSupply(recs, repoRoot, val('--supply', null), has('--history')); return; }
  if (has('--mcp-watch')) { cmdMcpWatch(recs); return; }
  if (has('--concurrency')) { cmdConcurrency(recs, repoRoot); return; }
  if (has('--m149-watch')) { cmdM149Watch(recs, repoRoot, has('--record')); return; }
  if (has('--m149')) { cmdM149(recs, repoRoot); return; }
  if (has('--diag')) {
    // ★digest は `idiom-recur.mjs --rescan` が作る。v3 でないと `ag`(toolUseId)が無い。
    const dp = val('--digest', path.join(repoRoot, '.cache', 'idiom-recur-digest.json'));
    if (!fs.existsSync(dp)) { console.error(`digest がない: ${dp}\n  まず \`node tools/idiom-recur.mjs --rescan\``); process.exit(2); }
    const dg = JSON.parse(fs.readFileSync(dp, 'utf8'));
    if ((dg.v || 1) < 3) { console.error(`digest が古い(v${dg.v || 1}、要 v3 以上)。\`node tools/idiom-recur.mjs --rescan\` で作り直すこと。`); process.exit(2); }
    cmdDiag(recs, dg);
    return;
  }
  if (has('--solo')) {
    // ★★M159: M139 の「他と混ざらない命令だけ」の**規則をコードで固定して**数え直す。
    const t = Date.now();
    const sc = await scanSoloCalls({ cache: path.join(repoRoot, '.cache', 'solo-timing.json') });
    const tal = tallySolo(sc.calls, [0, 1, 2, 3]);
    const LV = {
      0: '★L0 素(heredoc を剥がない)',
      1: '★★L1 正典 SOLO-1(= M153 の定義)',
      2: '★L2 厳格(他の tools/ も許さない)',
      3: '★L3 素朴(命令が 1 本だけ。`&&` `;` で繋がっていない)',
    };
    console.log('== --solo —— 「他と混ざらない命令だけ」(★規則はコードに焼いてある。M159) ==');
    console.log(`  ログ ${sc.nFiles} 本 / 親の Bash 呼び出し ${sc.bashUses} 件 / sidechain 除外 ${sc.sidechain} / 対にならず ${sc.unpaired}`);
    console.log(`  走査 ${((Date.now() - t) / 1000).toFixed(1)}s${sc.cached ? '(cache)' : ''}`);
    if (sc.calls.length) {
      const d0 = new Date(sc.calls[0].t0), d1 = new Date(sc.calls[sc.calls.length - 1].t0);
      console.log(`  期間 ${d0.toISOString().slice(0, 10)} 〜 ${d1.toISOString().slice(0, 10)}`
        + `(${((d1 - d0) / 86400000).toFixed(1)} 日)`);
    }
    for (const lv of [0, 1, 2, 3]) {
      const g = tal[lv];
      console.log('');
      console.log(`  ${LV[lv]}   —— 合計 ${g.n} 件 / ★${g.hours.toFixed(2)} 時間`);
      console.log('    道具              件数    合計h    中央値s');
      for (const r of g.rows) {
        console.log(`    ${r.id.padEnd(16)} ${String(r.n).padStart(5)} ${r.hours.toFixed(2).padStart(8)}`
          + `   ${(r.med === null ? '—' : r.med.toFixed(1)).padStart(7)}`);
      }
    }
    console.log('');
    console.log('  ★過去の公表値との比較(★同じ母集団ではないので「訂正」ではなく「定義の差」):');
    console.log(`     M139(メタ第 29 回。規則は台帳に無い)  合計 2.70 時間`);
    console.log(`     M153(メタ第 31 回。剥いだ・私の定義)  合計 4.52 時間`);
    console.log(`     ★M159 L1(正典。この道具)             合計 ${tal[1].hours.toFixed(2)} 時間`);
    return;
  }
  if (has('--main')) {
    // ★本体セッションの費用(M125)。★族の割り当てに digest を使う(無ければ族の表だけ空になる)。
    const dp = val('--digest', path.join(repoRoot, '.cache', 'idiom-recur-digest.json'));
    let dg = null;
    if (fs.existsSync(dp)) {
      dg = JSON.parse(fs.readFileSync(dp, 'utf8'));
      if ((dg.v || 1) < 3) { console.error(`digest が古い(v${dg.v || 1}、要 v3 以上)。\`node tools/idiom-recur.mjs --rescan\``); process.exit(2); }
    } else console.error(`※ digest が無い(${dp}) —— 族の表は空になる。\`node tools/idiom-recur.mjs --rescan\` で作れる。`);
    await cmdMain(recs, dg, {
      cache: path.join(repoRoot, '.cache', 'main-timing.json'),
      limit: Number(val('--limit', '15')),
    });
    return;
  }
  if (has('--cost')) {
    // ★worktree の decisions-pending.md は本体より古い(M99)。★--costs で本体の版を指せる。
    const cp = val('--costs', path.join(repoRoot, 'ResearchPaper', 'decisions-pending.md'));
    cmdCost(recs, [cp]);
    return;
  }

  // 既定: 一覧
  const lim = Number(val('--limit', '30'));
  const show = lim > 0 ? recs.slice(-lim) : recs;
  console.log('  完了時刻(UTC)          分   tool  tokens  type            持ち場');
  for (const r of show) {
    console.log(`  ${padr((r.ts ?? '').replace('T', ' ').slice(0, 19), 20)} ${pad(min1(r.durationMs), 6)} ${pad(r.toolUses, 5)} ${pad(r.tokens, 8)}  ${padr(r.agentType.slice(0, 14), 15)} ${String(r.description).slice(0, 44)}`);
  }
  if (lim > 0 && recs.length > lim) console.log(`  … 残り ${recs.length - lim} 件(--limit 0 で全部)`);
}

if (import.meta.url.endsWith(process.argv[1].replace(/\\/g, '/')) || process.argv[1].endsWith('agent-timing.mjs')) {
  main();
}
